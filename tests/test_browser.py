from pathlib import Path
import curses
import tempfile
import unittest
from unittest.mock import patch

import pysam

from tigv.core import (Tracks, View, bigwig_bins, coverage_and_reads, layout_reads,
                       parse_region, read_cells, read_connection_cells, read_fastq)
from tigv.app import Browser, gene_lanes
from tigv.gtf import Gene, GeneIndex

try:
    import pyBigWig
except ImportError:
    pyBigWig = None


class FakeScreen:
    def __init__(self):
        self.writes = []

    def getmaxyx(self):
        return 24, 80

    def erase(self):
        self.writes.clear()

    def addnstr(self, row, col, value, limit, attr):
        self.writes.append((row, col, value[:limit], attr))

    def refresh(self):
        pass

    def border(self):
        pass


class BrowserTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.fa = self.root / "ref.fa"
        self.fa.write_text(">chr1\n" + "A" * 200 + "\n")
        pysam.faidx(str(self.fa))
        self.gtf = self.root / "genes.gtf"
        self.gtf.write_text(
            'chr1\tx\tgene\t11\t80\t.\t+\t.\tgene_id "g1"; gene_name "ALPHA";\n'
            'chr1\tx\texon\t11\t20\t.\t+\t.\tgene_id "g1"; gene_name "ALPHA";\n'
            'chr1\tx\texon\t51\t60\t.\t+\t.\tgene_id "g1"; gene_name "ALPHA";\n'
        )
        self.bam = self.root / "reads.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 200}]}
        with pysam.AlignmentFile(str(self.bam), "wb", header=header) as output:
            read = pysam.AlignedSegment()
            read.query_name = "one"
            read.query_sequence = "AATAAA"
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 10
            read.mapping_quality = 60
            read.cigartuples = [(0, 3), (1, 1), (0, 1), (2, 1), (0, 1)]
            read.query_qualities = pysam.qualitystring_to_array("IIIIII")
            read.set_tag("MD", "2A1^A1")
            output.write(read)
            for start in (30, 50):
                mate = pysam.AlignedSegment()
                mate.query_name = "pair"
                mate.query_sequence = "AAAAAA"
                mate.flag = 0
                mate.reference_id = 0
                mate.reference_start = start
                mate.mapping_quality = 60
                mate.cigarstring = "6M"
                mate.query_qualities = pysam.qualitystring_to_array("IIIIII")
                output.write(mate)
            other = pysam.AlignedSegment()
            other.query_name = "one"
            other.query_sequence = "AAAAAA"
            other.flag = 2048  # supplementary alignment of the same read
            other.reference_id = 0
            other.reference_start = 100
            other.mapping_quality = 60
            other.cigarstring = "6M"
            other.query_qualities = pysam.qualitystring_to_array("IIIIII")
            output.write(other)
        pysam.index(str(self.bam))

    def tearDown(self):
        self.tmp.cleanup()

    def test_gene_search_and_exon_intron_structure(self):
        index = GeneIndex.from_file(self.gtf)
        gene = index.search("alpha")[0]
        self.assertEqual((gene.start, gene.end), (10, 80))
        self.assertEqual(gene.merged_exons(), [(10, 20), (50, 60)])
        self.assertEqual(index.overlap("chr1", 45, 55), [gene])
        self.assertEqual(index.search("g1"), [gene])

    def test_coordinate_parsing(self):
        self.assertEqual(parse_region("chr1:1,001-2,000"), ("chr1", 1000, 2000))
        self.assertEqual(parse_region("chr1:15"), ("chr1", 14, 15))
        with self.assertRaises(ValueError):
            parse_region("chr1:20-10")

    def test_browser_gene_and_region_navigation(self):
        tracks = Tracks()
        try:
            tracks.set_fasta(str(self.fa))
            tracks.set_gtf(str(self.gtf))
            browser = Browser(tracks, View("chr1", 0, 80))
            browser.search_gene("ALPHA", 40)
            self.assertEqual(browser.view.chrom, "chr1")
            self.assertLessEqual(browser.view.start, 10)
            self.assertGreaterEqual(browser.view.end, 80)
            browser.goto("1:100-120", 40)
            self.assertLessEqual(browser.view.start, 99)
            self.assertGreaterEqual(browser.view.end, 120)
        finally:
            tracks.close()

    def test_slash_search_jumps_to_locus_and_keeps_gene_and_read_search(self):
        tracks = Tracks()
        try:
            tracks.set_fasta(str(self.fa))
            tracks.set_gtf(str(self.gtf))
            tracks.add_bam(str(self.bam))
            browser = Browser(tracks, View("chr1", 0, 80))
            screen = FakeScreen()
            browser.prompt = lambda _screen, _prefix: "ALPHA"
            browser.key(screen, "/")
            self.assertEqual(browser.search_mode, "gene")

            browser.prompt = lambda _screen, _prefix: "1:100-120"
            browser.key(screen, "/")
            self.assertEqual(browser.view.chrom, "chr1")
            self.assertLessEqual(browser.view.start, 99)
            self.assertGreaterEqual(browser.view.end, 120)
            self.assertIsNone(browser.search_mode)
            browser.key(screen, "n")
            self.assertEqual(browser.status, "No previous gene or read search")

            browser.prompt = lambda _screen, _prefix: "chr1:15"
            browser.key(screen, "/")
            self.assertLessEqual(browser.view.start, 14)
            self.assertGreater(browser.view.end, 14)
            browser.prompt = lambda _screen, _prefix: "chr1"
            browser.key(screen, "/")
            self.assertEqual(browser.view.start, 0)

            browser.prompt = lambda _screen, _prefix: "@one"
            browser.key(screen, "/")
            self.assertEqual(browser.search_mode, "read")
            self.assertEqual(browser.selected_read_name, "one")
            browser.prompt = lambda _screen, _prefix: "chr1:20-10"
            browser.key(screen, "/")
            self.assertIn("Coordinates must be positive", browser.status)
        finally:
            tracks.close()

    def test_bam_mismatch_indel_and_coverage(self):
        tracks = Tracks()
        try:
            tracks.set_fasta(str(self.fa))
            tracks.add_bam(str(self.bam))
            view = View("chr1", 10, 20)
            coverage, reads, truncated = coverage_and_reads(tracks.bams[0][1], view, 20)
            self.assertFalse(truncated)
            self.assertEqual(len(reads), 1)
            self.assertGreater(max(coverage), 0)
            cells = read_cells(reads[0], view, 20, tracks.reference("chr1", 10, 30))
            self.assertEqual(cells[2], ("T", True))
            self.assertEqual(cells[3], ("+", True))
            self.assertEqual(cells[4], ("-", True))
            self.assertEqual(read_cells(reads[0], view, 20, None)[2], ("T", True))
            _, forced, limited = coverage_and_reads(tracks.bams[0][1], view, 20, max_reads=0,
                                                    highlighted_name="one")
            self.assertTrue(limited)
            self.assertEqual(len(forced), 1)
        finally:
            tracks.close()

    def test_indels_require_base_resolution_and_complete_visibility(self):
        header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 200}]})

        def make_read(cigar, sequence):
            read = pysam.AlignedSegment(header)
            read.query_name = "indel"
            read.query_sequence = sequence
            read.reference_id = 0
            read.reference_start = 10
            read.cigarstring = cigar
            return read

        deletion = make_read("2M3D2M", "AAAA")
        wide = View("chr1", 10, 20)
        coarse = read_cells(deletion, wide, 10, "A" * wide.span)
        self.assertIn("#", [glyph for glyph, _ in coarse.values()])
        self.assertNotIn("-", [glyph for glyph, _ in coarse.values()])
        self.assertNotIn("·", read_connection_cells([deletion], wide, 10).values())

        partial = View("chr1", 13, 5)
        self.assertNotIn("-", [glyph for glyph, _ in
                               read_cells(deletion, partial, 5, "A" * partial.span).values()])
        complete = View("chr1", 10, 7)
        complete_cells = read_cells(deletion, complete, 7, "A" * complete.span)
        self.assertEqual([complete_cells[col] for col in (2, 3, 4)], [("-", True)] * 3)

        insertion = make_read("2M6I2M", "A" * 10)
        short = View("chr1", 10, 4)
        self.assertNotIn("+", [glyph for glyph, _ in
                               read_cells(insertion, short, 4, "A" * short.span).values()])
        boundary = View("chr1", 10, 2)
        self.assertNotIn("+", [glyph for glyph, _ in
                               read_cells(insertion, boundary, 2, "A" * boundary.span).values()])
        full = View("chr1", 10, 10)
        self.assertEqual(read_cells(insertion, full, 10, "A" * full.span)[2], ("+", True))

    def test_coarse_mismatch_uses_star_and_match_uses_hash(self):
        tracks = Tracks()
        try:
            tracks.set_fasta(str(self.fa))
            tracks.add_bam(str(self.bam))
            view = View("chr1", 10, 20)
            _, reads, _ = coverage_and_reads(tracks.bams[0][1], view, 10)
            cells = read_cells(reads[0], view, 10, tracks.reference("chr1", 10, 30))
            glyphs = {glyph for glyph, _ in cells.values()}
            self.assertIn("#", glyphs)
            self.assertIn("*", glyphs)
            self.assertNotIn("+", glyphs)
            self.assertNotIn("-", glyphs)
        finally:
            tracks.close()

    def test_fastq_reconstruction_uses_original_read_orientation(self):
        read = pysam.AlignedSegment()
        read.query_name = "reverse"
        read.flag = 16
        read.query_sequence = "ACTT"
        read.query_qualities = pysam.qualitystring_to_array("abcd")
        self.assertEqual(read_fastq(read), ("@reverse\nAAGT\n+\ndcba", None))
        read.query_qualities = None
        fastq, reason = read_fastq(read)
        self.assertIsNone(fastq)
        self.assertIn("qualities", reason)

    def test_mouse_click_opens_read_detail_window_with_fastq(self):
        tracks = Tracks()
        try:
            tracks.set_fasta(str(self.fa))
            tracks.add_bam(str(self.bam))
            browser = Browser(tracks, View("chr1", 0, 80))
            screen = FakeScreen()
            browser.draw(screen)
            y, x = next((y, x) for (y, x), selection in browser.read_click_map.items()
                        if x >= 13 and selection.read.query_name == "one")
            with patch("tigv.app.curses.getmouse", return_value=(0, x, y, 0, curses.BUTTON1_CLICKED)):
                browser.key(screen, curses.KEY_MOUSE)
            self.assertEqual(browser.read_detail.read.query_name, "one")
            self.assertIn(("@one", False), browser.detail_lines)
            self.assertIn(("AATAAA", True), browser.detail_lines)
            self.assertIn(("IIIIII", True), browser.detail_lines)
            popup = FakeScreen()
            with patch("tigv.app.curses.newwin", return_value=popup) as newwin:
                browser.draw(screen)
            newwin.assert_called_once()
            self.assertTrue(any("Read details" in value for _, _, value, _ in popup.writes))
            self.assertTrue(any("FASTQ" in value for _, _, value, _ in popup.writes))
            browser.key(screen, "q")
            self.assertIsNone(browser.read_detail)
            self.assertTrue(browser.key(screen, "r"))
        finally:
            tracks.close()

    def test_read_name_search_jump_highlight_and_track_focus(self):
        tracks = Tracks()
        try:
            tracks.set_fasta(str(self.fa))
            tracks.set_gtf(str(self.gtf))
            tracks.add_bam(str(self.bam))
            hits = tracks.find_read("one")
            self.assertEqual([hit.start for hit in hits], [10, 100])
            self.assertIs(hits, tracks.find_read("one"))
            browser = Browser(tracks, View("chr1", 0, 80))
            browser.search_read("one", 66)
            self.assertEqual(browser.focus_index, 2)  # Reference, Genes, BAM
            self.assertEqual(browser.selected_read_name, "one")
            self.assertLessEqual(browser.view.start, 10)
            self.assertGreater(browser.view.end, 10)
            screen = FakeScreen()
            browser.draw(screen)
            self.assertFalse(browser.pending_read_reveal)
            self.assertTrue(any(value.startswith("*one") for _, _, value, _ in screen.writes))
            browser.key(screen, "n")
            self.assertEqual(browser.view.chrom, "chr1")
            self.assertLessEqual(browser.view.start, 100)
            self.assertGreater(browser.view.end, 100)
            browser.key(screen, "N")
            self.assertEqual(browser.read_hit_index, 0)
            browser.key(screen, "K")
            self.assertEqual(browser.focus_index, 1)
            browser.key(screen, "K")
            self.assertEqual(browser.focus_index, 0)
            browser.key(screen, "J")
            self.assertEqual(browser.focus_index, 1)
            browser.key(screen, "e")
            self.assertTrue(browser.gene_expanded)
            browser.command("collapse", 66)
            self.assertFalse(browser.gene_expanded)
            browser.command("expand", 66)
            self.assertTrue(browser.gene_expanded)
            browser.command("clear", 66)
            self.assertIsNone(browser.selected_read_name)
            browser.command("read one", 66)
            self.assertEqual(browser.selected_read_name, "one")
            browser._search("@one", 66)
            self.assertEqual(browser.search_mode, "read")
        finally:
            tracks.close()

    def test_expanded_gene_rows_never_overlap(self):
        genes = [Gene("chr1", 10, 80, "g1", "A"), Gene("chr1", 15, 50, "g2", "B"),
                 Gene("chr1", 90, 120, "g3", "C")]
        view = View("chr1", 0, 140)
        packed = gene_lanes(genes, view, 80, False)
        expanded = gene_lanes(genes, view, 80, True)
        self.assertEqual(len(packed), 2)
        self.assertEqual(len(expanded), 3)
        self.assertTrue(all(len(row) == 1 for row in expanded))
        tracks = Tracks()
        tracks.genes = GeneIndex(genes)
        browser = Browser(tracks, view)
        browser.gene_expanded = True
        screen = FakeScreen()
        browser.draw(screen)
        rows = {row for row, col, value, _ in screen.writes
                if col == 0 and value.strip() in {"A", "B", "C"}}
        self.assertEqual(len(rows), 3)

    def test_same_read_pieces_share_a_row_and_are_joined_by_dots(self):
        tracks = Tracks()
        try:
            tracks.set_fasta(str(self.fa))
            tracks.add_bam(str(self.bam))
            browser = Browser(tracks, View("chr1", 0, 80))
            browser.search_read("pair", 66)
            screen = FakeScreen()
            browser.draw(screen)
            self.assertTrue(any(value == "·" for _, _, value, _ in screen.writes))
            _, reads, _ = coverage_and_reads(tracks.bams[0][1], View("chr1", 0, 80), 80)
            rows = layout_reads(reads, View("chr1", 0, 80), 80)
            self.assertTrue(any(sum(read.query_name == "pair" for read in row) == 2 for row in rows))
        finally:
            tracks.close()

    def test_supplementary_alignment_cross_gene_renders_brackets_on_one_row(self):
        cross_gtf = self.root / "cross.gtf"
        cross_gtf.write_text(
            'chr1\tx\tgene\t11\t20\t.\t+\t.\tgene_id "left"; gene_name "LEFT";\n'
            'chr1\tx\texon\t11\t20\t.\t+\t.\tgene_id "left"; gene_name "LEFT";\n'
            'chr1\tx\tgene\t101\t110\t.\t+\t.\tgene_id "right"; gene_name "RIGHT";\n'
            'chr1\tx\texon\t101\t110\t.\t+\t.\tgene_id "right"; gene_name "RIGHT";\n'
        )
        tracks = Tracks()
        try:
            tracks.set_gtf(str(cross_gtf))
            tracks.add_bam(str(self.bam))
            view = View("chr1", 0, 120)
            _, reads, _ = coverage_and_reads(tracks.bams[0][1], view, 66)
            self.assertEqual(sum(read.query_name == "one" for read in reads), 2)
            rows = layout_reads(reads, view, 66, tracks.genes)
            self.assertTrue(any(sum(read.query_name == "one" for read in row) == 2 for row in rows))
            browser = Browser(tracks, view)
            screen = FakeScreen()
            browser.draw(screen)
            left_rows = {row for row, _, value, _ in screen.writes if value == "["}
            right_rows = {row for row, _, value, _ in screen.writes if value == "]"}
            self.assertTrue(left_rows & right_rows)
        finally:
            tracks.close()

    def test_exons_join_with_dots_and_cross_gene_chrna_uses_brackets(self):
        header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 200_000}]})

        def make_read(name, start, cigar="10M", sequence="A" * 10):
            read = pysam.AlignedSegment(header)
            read.query_name = name
            read.query_sequence = sequence
            read.reference_id = 0
            read.reference_start = start
            read.cigarstring = cigar
            return read

        near_left = make_read("near", 10)
        near_right = make_read("near", 30)
        near = read_connection_cells([near_left, near_right], View("chr1", 0, 80), 80)
        self.assertEqual("".join(near[column] for column in range(20, 30)), "·" * 10)
        spliced = make_read("splice", 10, "10M20N10M", "A" * 20)
        self.assertIn("·", read_connection_cells([spliced], View("chr1", 0, 80), 80).values())

        far_left = make_read("far", 100)
        far_right = make_read("far", 150_000)
        same_gene = GeneIndex([Gene("chr1", 90, 150_010, "g1", "SAME",
                                    exons=[(100, 110), (150_000, 150_010)])])
        distant_same_gene = read_connection_cells([far_left, far_right],
                                                   View("chr1", 0, 200_000), 80, same_gene)
        self.assertIn("·", distant_same_gene.values())
        self.assertNotIn("[", distant_same_gene.values())
        self.assertNotIn("]", distant_same_gene.values())

        long_splice = make_read("long-splice", 100, "10M149890N10M", "A" * 20)
        splice_link = read_connection_cells([long_splice], View("chr1", 0, 200_000), 80, same_gene)
        self.assertIn("·", splice_link.values())
        self.assertNotIn("[", splice_link.values())

        two_genes = GeneIndex([
            Gene("chr1", 90, 111, "g1", "LEFT", exons=[(100, 110)]),
            Gene("chr1", 125, 150_010, "g2", "RIGHT", exons=[(130, 140), (150_000, 150_010)]),
        ])
        near_cross_gene = read_connection_cells([make_read("chimera", 100), make_read("chimera", 130)],
                                                View("chr1", 90, 60), 60, two_genes)
        drawing = "".join(near_cross_gene.get(column, " ") for column in range(60))
        self.assertIn("[", drawing)
        self.assertIn("]", drawing)
        self.assertNotIn("-", drawing)
        self.assertIn(" ", drawing[drawing.index("[") + 1 : drawing.index("]")])
        self.assertNotIn("·", near_cross_gene.values())

        far_cross_gene = read_connection_cells([far_left, far_right],
                                               View("chr1", 0, 200_000), 80, two_genes)
        self.assertIn("[", far_cross_gene.values())
        self.assertIn("]", far_cross_gene.values())

        lone = make_read("lone", 100)
        lone.set_tag("SA", "chr1,150001,+,10M,60,0;")
        open_end = read_connection_cells([lone], View("chr1", 0, 200), 80, two_genes)
        self.assertIn("[", open_end.values())
        self.assertNotIn("-", open_end.values())
        opposite = make_read("opposite", 150_000)
        opposite.set_tag("SA", "chr1,101,+,10M,60,0;")
        left_end = read_connection_cells([opposite], View("chr1", 149_900, 200), 80, two_genes)
        self.assertIn("]", left_end.values())
        self.assertNotIn("-", left_end.values())

        unknown = read_connection_cells([far_left, far_right], View("chr1", 0, 200_000), 80)
        self.assertIn("·", unknown.values())
        self.assertNotIn("[", unknown.values())

    @unittest.skipIf(pyBigWig is None, "pyBigWig unavailable")
    def test_bigwig_signal(self):
        path = self.root / "signal.bw"
        bw = pyBigWig.open(str(path), "w")
        bw.addHeader([("chr1", 200)])
        bw.addEntries(["chr1"], [10], ends=[20], values=[5.0])
        bw.close()
        tracks = Tracks()
        try:
            tracks.add_bigwig(str(path))
            values = bigwig_bins(tracks.bigwigs[0][1], View("chr1", 0, 40), 20)
            self.assertGreater(max(value for value in values if value is not None), 0)
        finally:
            tracks.close()


if __name__ == "__main__":
    unittest.main()
