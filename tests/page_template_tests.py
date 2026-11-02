"""Tests that pin the source text of the results server's page templates.

The templates are Jinja + JavaScript, and there is no Python harness that can execute the
JavaScript. These tests pin the template source text instead (rendering a Jinja macro where its
output is what matters), which is enough to catch a silent change to the behavior each test
describes.
"""

import os
import re
import unittest

import jinja2

CODE_DIRECTORY = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "trails")


def read_template(name):
    """Returns the text of a file in the trails/ code directory.

    Args:
        name: File name, e.g. "index_page_template.html" or "results_server.py".

    Returns:
        The file's contents as a string.
    """
    with open(os.path.join(CODE_DIRECTORY, name)) as handle:
        return handle.read()


def js_string_list(text, name):
    """Returns the string literals of a `const <name> = [...]` JavaScript array.

    Args:
        text: Template source text to search.
        name: Name of the const, e.g. "MOTIF_BINS".

    Returns:
        List of the single-quoted strings inside the array literal.
    """
    match = re.search(r"const %s = \[([^\]]*)\]" % name, text)
    assert match, "no `const %s = [...]` in the template" % name
    return re.findall(r"'([^']*)'", match.group(1))


def server_motif_bins():
    """Returns the motif-size bin names results_server.compute_sample_qc_data emits."""
    server = read_template("results_server.py")
    match = re.search(r"WHEN MotifSize IS NULL.*?END", server, re.DOTALL)
    assert match, "the motif-size CASE expression moved out of compute_sample_qc_data"
    return re.findall(r"(?:THEN|ELSE) '([^']*)'", match.group(0))


def js_function_body(text, signature):
    """Returns the source of a top-level JavaScript function, braces included.

    Args:
        text: Template source text to search.
        signature: Start of the function's declaration, e.g. "async function saveNote(".

    Returns:
        The text from the declaration through the closing brace in column zero.
    """
    start = text.index(signature)
    end = text.index("\n}\n", start)
    return text[start:end + 3]


def js_object_body(text, signature):
    """Returns the source of a top-level `const <name> = {...};` JavaScript object literal.

    Args:
        text: Template source text to search.
        signature: Start of the declaration, e.g. "const DEFAULT_STATE = {".

    Returns:
        The text from the declaration through the `};` in column zero that ends it.
    """
    start = text.index(signature)
    end = text.index("\n};\n", start)
    return text[start:end + 4]


def filter_control_macros():
    """Returns the macros of filter_control_macros.html, callable from Python.

    Returns:
        The template module, with gene_region_checkboxes(name_prefix) and
        population_filter_field(name_prefix, checkbox_attribute) as attributes.
    """
    return jinja2.Environment(loader=jinja2.FileSystemLoader(CODE_DIRECTORY)).get_template(
        "filter_control_macros.html").module


class TagFilterRestoreTests(unittest.TestCase):
    """The saved tag must be selected in a dropdown whose menu already exists.

    applyStateToUI selects the tag by value, and Semantic UI silently ignores a 'set selected'
    for a value the menu does not hold yet, so the Search and Export handlers (which re-read the
    dropdown through readStateFromUI) would drop the tag on the next request.
    """

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_startup_waits_for_both_dropdowns_before_restoring_url_state(self):
        self.assertIn("Promise.all([tagDropdownReady, sampleIdDropdownReady]).then(() => {\n"
                      "        // Read URL state and perform initial search after dropdowns are populated\n"
                      "        readStateFromUrl();", self.template)

    def test_the_dropdown_promises_are_kept_rather_than_fired_and_forgotten(self):
        for declaration in ("const tagDropdownReady = setUpTagFilterDropdown($('#tag-filter'));",
                            "const sampleIdDropdownReady = fetchJsonOrNull(`${API_BASE_URL}/schema`, 'Load schema')"):
            self.assertIn(declaration, self.template)
            self.assertLess(self.template.index(declaration), self.template.index("Promise.all([tagDropdownReady,"))
        # setUpTagFilterDropdown (header_template.html) hands back the fill's promise
        self.assertIn("    return refreshTagFilterDropdown($dropdown);\n}",
                      js_function_body(read_template("header_template.html"), "function setUpTagFilterDropdown("))

    def test_search_and_export_read_the_tag_off_the_dropdown(self):
        # This is what makes the ordering matter, so pin it: if these ever stopped reading the
        # dropdown the startup ordering above would no longer be load bearing.
        self.assertIn("currentState.tag = getDropdownValue($('#tag-filter'));", self.template)
        for handler in ("$('#search-button').on('click', function() {", "async function performExport(format) {"):
            self.assertIn("readStateFromUI();", self.template[self.template.index(handler):][:2000])


class SortIconTests(unittest.TestCase):
    """Every sortable column header needs the span the active-sort styling colors."""

    def test_all_sortable_headers_carry_a_sort_icon(self):
        template = read_template("index_page_template.html")
        headers = re.findall(r"<th[^>]*\bdata-sort=\"([^\"]+)\"[^>]*>(.*?)</th>", template)
        self.assertEqual({"region", "size", "count", "family_count", "variation_cluster"},
                         {field for field, _ in headers})
        for field, body in headers:
            self.assertIn('<span class="sort-icon">', body, field)

    def test_the_active_sort_class_is_only_styled_through_that_span(self):
        # Why the span is required: header_template.html colors .sort-icon, nothing else.
        self.assertIn("#results-table thead th.sorted-active .sort-icon", read_template("header_template.html"))


class ActiveSortClassTests(unittest.TestCase):
    """The sort indicator class must not claim a direction the server does not apply.

    A header click only adds or removes the field from the sort list; each field's direction is
    fixed in results_server.build_api_order_by's SORT_MAPPING, and Region sorts ascending, so one
    direction-free class marks every sorted column.
    """

    def test_the_class_the_indicator_applies_is_the_class_the_stylesheet_matches(self):
        indicator = js_function_body(read_template("index_page_template.html"),
                                     "function syncSortHeaderIndicator(")
        self.assertIn("removeClass('sorted-active')", indicator)
        self.assertIn("addClass('sorted-active')", indicator)
        styled = re.findall(r"#results-table thead th\.([\w-]+) \.sort-icon",
                            read_template("header_template.html"))
        self.assertEqual(["sorted-active"], styled)

    def test_the_server_sorts_the_region_column_ascending(self):
        # Region sorts ascending while the other columns sort descending, so a sorted-column
        # marker cannot name a direction.
        self.assertIn('"region": ("gene_region_rank", "ASC")', read_template("results_server.py"))


class DetailRowCellSpecificityTests(unittest.TestCase):
    """The detail-row cell rule must not reach the tables nested inside that cell.

    A descendant selector, `#results-table tbody tr.locus-detail-row td` (1 id, 1 class, 3 types),
    would out-specify `#results-table .sample-subtable td` (1,1,1), `.threshold-marker` (0,1,0) and
    `.invisible-table td` (0,1,1), giving every sample cell the outer cell's 40px indent, #fafafa
    background and 2px bottom border.
    """

    def setUp(self):
        self.header = read_template("header_template.html")

    def test_the_detail_cell_rule_uses_a_child_combinator(self):
        self.assertIn("#results-table tbody tr.locus-detail-row > td {", self.header)
        self.assertNotIn("#results-table tbody tr.locus-detail-row td {", self.header)

    def test_the_detail_row_has_exactly_one_direct_cell(self):
        # The child combinator is only equivalent for the outer cell because the row is built
        # with a single colspan td that the detail HTML is dropped into.
        template = read_template("index_page_template.html")
        self.assertIn('<tr class="locus-detail-row"><td colspan="${$row.children(\'td\').length}"', template)
        self.assertEqual(1, template.count('<tr class="locus-detail-row">'))

    def test_the_invisible_table_rule_outranks_the_main_table_cell_rule(self):
        # `.invisible-table td` alone loses to `#results-table tbody td`, which the Affected
        # Unsolved summary table sits inside.
        self.assertIn(".invisible-table td,\n#results-table .invisible-table td {", self.header)
        self.assertIn("#results-table tbody td {", self.header)

    def test_the_compact_sub_table_rules_win(self):
        self.assertIn("#results-table .sample-subtable td { padding: 2px 6px;", self.header)


class AnnotationWriteHelperTests(unittest.TestCase):
    """The results page and the swim plot share one set of annotation write helpers.

    They are defined once, in header_template.html, which both pages include, and every write
    failure surfaces the server's own explanation from the response body instead of a bare status.
    """

    PAGES_THAT_WRITE_ANNOTATIONS = ("index_page_template.html", "swim_plot_template.html")
    HELPER_SIGNATURES = ("async function saveNote(", "async function deleteNote(",
                         "async function addTag(", "async function removeTag(")

    def setUp(self):
        self.header = read_template("header_template.html")

    def test_both_pages_use_the_one_shared_copy(self):
        for name in self.PAGES_THAT_WRITE_ANNOTATIONS:
            template = read_template(name)
            self.assertIn('{% include "header_template.html" %}', template, name)
            for signature in self.HELPER_SIGNATURES:
                self.assertIn(signature, self.header)
                self.assertNotIn(signature, template, name)
            for call in ("openNoteEditor(", "addTag(", "removeTag("):
                self.assertIn(call, template, name)
        self.assertIn("saveNote(currentLocusId, noteText)",
                      js_function_body(read_template("note_editor.html"), "async function openNoteEditor("))

    def test_the_write_helpers_read_the_response_body(self):
        for signature in self.HELPER_SIGNATURES:
            self.assertIn("throw await errorFromResponse(response);", js_function_body(self.header, signature), signature)

    def test_no_template_names_a_server_symbol_that_does_not_exist(self):
        # Every `results_server.<name>` a comment cites has to be findable in the server.
        server = read_template("results_server.py")
        for name in sorted(os.listdir(CODE_DIRECTORY)):
            if not name.endswith(".html"):
                continue
            for symbol in set(re.findall(r"results_server\.(?!py\b)([a-zA-Z_][a-zA-Z_0-9]*)", read_template(name))):
                self.assertTrue(re.search(r"^(def %s\(|%s = )" % (symbol, symbol), server, re.MULTILINE),
                                "%s cites results_server.%s, which does not exist" % (name, symbol))

    def test_the_read_only_server_hook_answers_403_with_a_detail_field(self):
        server = read_template("results_server.py")
        self.assertIn("def restrict_writes():", server)
        hook = server[server.index("def restrict_writes():"):]
        hook = hook[:hook.index("\n@app.")]
        self.assertIn('"detail": "This server is bound to a non-loopback address', hook)
        self.assertIn("403", hook)

    def test_the_server_sends_the_fields_the_helpers_read(self):
        server = read_template("results_server.py")
        self.assertIn('"error": "annotation changes are limited to this machine"', server)
        self.assertIn('{"error": "note_text is required and cannot be empty"}', server)


class OutlierWarningTests(unittest.TestCase):
    """A sample's warning icon combines its Mendelian-violation and its outlier warning.

    The QC2 page shows only the Mendelian part, so it leaves the outlier map out of the page.
    """

    def test_the_pages_that_show_outlier_warnings_define_the_map(self):
        for name in ("index_page_template.html", "sample_qc_template.html", "swim_plot_template.html"):
            self.assertIn("const OUTLIER_WARNINGS = {{ outlier_warnings | tojson }};", read_template(name), name)
        # The index, swim and sample QC routes render it
        self.assertEqual(3, read_template("results_server.py").count(
            'outlier_warnings=app.config.get("OUTLIER_WARNINGS", {}),'))

    def test_qc2_page_asks_for_the_mendelian_part_only(self):
        template = read_template("qc2_template.html")
        self.assertEqual(["item.sample_id, { mendelianOnly: true }"],
                         re.findall(r"renderMendelianWarning\(([^)]+)\)", template))
        self.assertNotIn("outlier_warnings", template)

    def test_the_popup_shows_the_motif_bin_labels_the_server_sends(self):
        # The server keys motif_rates by MENDELIAN_MOTIF_SIZE_BINS's labels, so the popup keeps no
        # label table of its own
        self.assertIn('warning["motif_rates"][label] = {', read_template("results_server.py"))
        self.assertIn("for (const [label, data] of Object.entries(warning.motif_rates)) {",
                      read_template("qc_shared_js.html"))

    def test_the_helper_reads_the_outlier_map_only_without_mendelian_only(self):
        # Leaving the map out of the QC2 page is only safe because the helper does not read it then
        helper = js_function_body(read_template("qc_shared_js.html"), "function renderMendelianWarning(")
        self.assertIn("function renderMendelianWarning(sampleId, { mendelianOnly = false } = {}) {", helper)
        self.assertIn("const outlier = mendelianOnly ? undefined : OUTLIER_WARNINGS[sampleId];", helper)


class MendelianWarningPopupTests(unittest.TestCase):
    """The warning popup highlights the rates the server flags a sample for, and only those.

    load_mendelian_warnings flags a sample when its autosome or chrX rate is above
    MENDELIAN_VIOLATION_RATE_WARNING_THRESHOLD; its chrY and motif-size rates are shown for context.
    """

    def setUp(self):
        self.shared_js = read_template("qc_shared_js.html")
        self.server = read_template("results_server.py")

    def test_the_highlight_threshold_is_the_servers(self):
        self.assertIn("MENDELIAN_VIOLATION_RATE_WARNING_THRESHOLD = 0.10", self.server)
        self.assertIn("const formatRate = rate => rate > 0.10 ? ", self.shared_js)

    def test_only_the_rows_the_server_flags_on_are_highlighted(self):
        for category in ("autosome", "chrX"):
            self.assertIn('rates["%s"] > MENDELIAN_VIOLATION_RATE_WARNING_THRESHOLD' % category, self.server)
            self.assertIn("${formatRate(warning.%s_rate)}" % category, self.shared_js)
        self.assertEqual(2, self.shared_js.count("${formatRate("))
        self.assertIn("${formatPercent(warning.chrY_rate)}", self.shared_js)
        self.assertIn("${formatPercent(data.rate)}", self.shared_js)


class ClipboardFallbackTests(unittest.TestCase):
    """The copy button must survive an unavailable or refused clipboard API."""

    def setUp(self):
        self.handler = read_template("index_page_template.html")
        start = self.handler.index("$(document).on('click', '.copy-ids-btn'")
        self.handler = self.handler[start:self.handler.index("\n});\n", start)]

    def test_availability_is_checked_before_the_write(self):
        self.assertIn("if (!navigator.clipboard || !navigator.clipboard.writeText) {", self.handler)
        self.assertLess(self.handler.index("!navigator.clipboard"),
                        self.handler.index("navigator.clipboard.writeText(text)"))

    def test_a_rejected_write_falls_back_instead_of_going_unhandled(self):
        self.assertIn(".catch(() => copyByHand());", self.handler)

    def test_the_fallback_shows_the_ids_the_user_asked_for(self):
        self.assertIn("const copyByHand = () => window.prompt('Copy the affected sample IDs:', text);",
                      self.handler)


class ExpandHandlerScopeTests(unittest.TestCase):
    """The expand/collapse handler must not fire for the table's header cell."""

    def test_handler_is_bound_to_table_cells_only(self):
        template = read_template("index_page_template.html")
        self.assertIn("$(document).on('click', 'td.col-expand'", template)
        self.assertNotIn("$(document).on('click', '.col-expand'", template)

    def test_header_cell_carries_the_class(self):
        # The th keeps the class for column sizing, which is why the handler has to be scoped.
        self.assertIn('<th class="col-expand">', read_template("index_page_template.html"))


class TagAttributeReadTests(unittest.TestCase):
    """A tag named "true" or "1" must survive being read back off the badge.

    jQuery's .data() converts such attribute values to a boolean or a number, which then fails
    the strict comparison against the active tag filter, so the page never reloads.
    """

    def test_results_page_reads_the_tag_attribute_as_a_string(self):
        template = read_template("index_page_template.html")
        self.assertIn("const tag = $(this).attr('data-tag');", template)
        self.assertNotIn(".data('tag')", template)

    def test_swim_plot_reads_the_tag_attribute_as_a_string(self):
        template = read_template("swim_plot_template.html")
        self.assertIn("var tag = $(this).attr('data-tag');", template)
        self.assertNotIn(".data('tag')", template)


class SharedJsHeaderTests(unittest.TestCase):
    """qc_shared_js.html's header comment must match what the file defines and who uses it."""

    def setUp(self):
        self.shared_js = read_template("qc_shared_js.html")
        self.header = self.shared_js.split("=========================================================================")[1]

    def test_header_names_every_including_template(self):
        for name in sorted(os.listdir(CODE_DIRECTORY)):
            if not name.endswith(".html") or name == "qc_shared_js.html":
                continue
            if '{% include "qc_shared_js.html" %}' in read_template(name):
                self.assertIn(name, self.header)

    def test_header_names_every_helper_the_file_defines(self):
        helpers = re.findall(r"^function (\w+)\(", self.shared_js, re.MULTILINE)
        self.assertIn("isAffectedStatus", helpers)
        self.assertIn("setupDropdownMenu", helpers)
        for helper in helpers:
            self.assertIn(helper, self.header)

    def test_header_says_only_mendelian_warnings_is_always_required(self):
        # QC2 defines no OUTLIER_WARNINGS (OutlierWarningTests), so the header must not ask for it
        self.assertIn("Requires Chart.js + jQuery (loaded by header_template.html) and the page-level "
                      "MENDELIAN_WARNINGS\n// const.", self.header)
        self.assertIn("mendelianOnly: true (QC2 does, and does not define it)", self.header)


class UnknownMotifBinTests(unittest.TestCase):
    """The 'Unknown' motif-size bin the server emits must reach the charts and the popup."""

    def test_the_server_emits_an_unknown_bin(self):
        self.assertIn("Unknown", server_motif_bins())

    def test_sample_qc_charts_cover_every_bin_the_server_emits(self):
        template = read_template("sample_qc_template.html")
        bins = set(js_string_list(template, "MOTIF_BINS"))
        bins.add(re.search(r"const UNKNOWN_MOTIF_BIN = '([^']*)'", template).group(1))
        self.assertEqual(bins, set(server_motif_bins()))

    def test_unknown_bin_is_rendered_only_when_the_data_has_it(self):
        template = read_template("sample_qc_template.html")
        self.assertIn("if (data.some(item => item.bin === UNKNOWN_MOTIF_BIN)) bins.push(UNKNOWN_MOTIF_BIN);",
                      template)
        # The grouping and the chart loop must both run over the extended list, not MOTIF_BINS.
        self.assertIn("const binData = groupRowsByBin(data, bins);", template)
        self.assertIn("bins.forEach((bin, idx) => {", template)

    def test_warning_popup_bin_order_covers_every_bin_the_server_emits(self):
        self.assertEqual(set(js_string_list(read_template("qc_shared_js.html"), "binOrder")),
                         set(server_motif_bins()))


class QCPageSharedCodeTests(unittest.TestCase):
    """The Sample QC and QC2 pages share their stylesheet and their per-bin chart code.

    The .qc-* rules are in header_template.html and the per-bin card and box-plot chart builders
    in qc_shared_js.html, which both pages include, so each page keeps only its own chart data and
    options.
    """

    QC_PAGES = ("sample_qc_template.html", "qc2_template.html")

    def test_the_qc_styles_are_in_the_shared_header(self):
        self.assertIn(".qc-chart-canvas-wrapper {", read_template("header_template.html"))
        for name in self.QC_PAGES:
            self.assertNotIn("<style>", read_template(name), name)

    def test_every_bin_chart_goes_through_the_shared_builders(self):
        for name in self.QC_PAGES:
            template = read_template(name)
            self.assertIn('{% include "qc_shared_js.html" %}', template, name)
            self.assertIn("appendQCBinCard(container, idx, bin, ", template, name)
            self.assertIn("createSingleBinBoxPlotChart(canvasId, ", template, name)
            self.assertNotIn("new Chart(", template, name)
            self.assertNotIn("document.createElement(", template, name)

    def test_qc2_builds_every_section_with_one_function(self):
        # The per-trio and per-center sections differ only in TRIO_SECTION and CENTER_SECTION, so
        # the zero-denominator exclusion and the shared y-axis rounding exist once.
        template = read_template("qc2_template.html")
        self.assertEqual(1, template.count("data.filter(item => item.total > 0)"))
        self.assertEqual(1, template.count("Math.ceil(maxFraction / 0.05) * 0.05 || 0.05"))
        draw = template[template.index("drawQCPageOnLoad('qc2_data', "):]
        self.assertEqual(2, draw.count("_BINS, TRIO_SECTION);"))
        self.assertEqual(2, draw.count("_BINS, CENTER_SECTION);"))

    def test_both_pages_load_and_group_their_data_through_the_shared_helpers(self):
        # drawQCPageOnLoad holds the fetch, the loading indicator and the error message, and
        # groupRowsByBin the grouping of a section's rows, once for both pages
        for name, call in (("sample_qc_template.html", "drawQCPageOnLoad('sample_qc_data', 'sample QC data', "
                                                       "['section-rank1', 'section-top10'], data => {"),
                           ("qc2_template.html", "drawQCPageOnLoad('qc2_data', 'QC2 data', ['section-motif', ")):
            template = read_template(name)
            self.assertIn(call, template, name)
            self.assertIn("const binData = groupRowsByBin(", template, name)
            self.assertNotIn("fetch(", template, name)
            self.assertNotIn("DOMContentLoaded", template, name)
        draw = js_function_body(read_template("qc_shared_js.html"), "function drawQCPageOnLoad(")
        self.assertIn("document.addEventListener('DOMContentLoaded', async function() {", draw)
        self.assertIn("document.getElementById('qc-error').textContent = `Failed to load ${dataName}: ` + err.message;", draw)

    def test_the_center_section_titles_name_the_server_trio_minimum(self):
        # The cutoff is the server's MIN_TRIOS_PER_CENTER, which the qc2() route passes in
        template = read_template("qc2_template.html")
        self.assertEqual(2, template.count("centers with ≥{{ min_trios_per_center }} trios</div>"))
        self.assertIn("min_trios_per_center=MIN_TRIOS_PER_CENTER,", read_template("results_server.py"))

    def test_the_shared_builder_hands_the_box_plot_plugin_its_statistics(self):
        # Both take the values sorted once
        builder = js_function_body(read_template("qc_shared_js.html"), "function createSingleBinBoxPlotChart(")
        self.assertEqual(1, builder.count(".sort("))
        self.assertIn("chart.config._boxStats = computeBoxStats(sortedValues);", builder)
        self.assertIn("p95: computePercentile(sortedValues, 95),", builder)
        self.assertIn("p99: computePercentile(sortedValues, 99)", builder)

    def test_the_median_and_quartiles_come_from_one_median_helper(self):
        stats = js_function_body(read_template("qc_shared_js.html"), "function computeBoxStats(")
        self.assertIn("const median = medianOfSorted(sorted);", stats)
        self.assertIn("const q1 = medianOfSorted(sorted.slice(0, Math.floor(n / 2))) ?? median;", stats)
        self.assertIn("const q3 = medianOfSorted(sorted.slice(Math.ceil(n / 2))) ?? median;", stats)
        self.assertNotIn(".sort(", stats)

    def test_the_percentile_labels_have_room_inside_the_canvas(self):
        # The plugin draws the '95%' / '99%' labels (about 18px wide) 2px right of the chart area,
        # so the builder has to leave at least that much room there.
        template = read_template("qc_shared_js.html")
        self.assertIn("ctx.fillText(config.label, chartArea.right + 2, yPixel + 3);", template)
        builder = js_function_body(template, "function createSingleBinBoxPlotChart(")
        self.assertIn("layout: { padding: { right: 22 } },", builder)


class PopulationDataFilterGatingTests(unittest.TestCase):
    """include_loci_without_population_data only counts when require_above_population is set."""

    def test_the_results_query_gates_the_flag(self):
        # buildApiUrl, whose query both /loci and the locus export use, sends the checkbox flags
        # from CHECKBOX_FILTERS, and this one only under a require_above_population test.
        self.assertIn("    for (const name of CHECKBOX_FILTERS) {\n"
                      "        if (currentState[name]\n"
                      "                && (name !== 'include_loci_without_population_data' || currentState.require_above_population)) {\n"
                      "            params.set(name, 'true');\n"
                      "        }\n"
                      "    }", js_function_body(read_template("index_page_template.html"), "function buildApiUrl("))

    def test_the_filter_chip_is_suppressed_with_no_population_threshold(self):
        template = read_template("index_page_template.html")
        self.assertIn("if (key === 'include_loci_without_population_data' && !currentState.require_above_population) continue;",
                      template)
        # population_metric is never suppressed: the Δpop column and allele colours use it even with no
        # population filter, so a non-default metric always changes what the table shows
        self.assertNotIn("if (key === 'population_metric'", js_function_body(template, "function updateActiveFiltersBar("))


class LocusExportUrlTests(unittest.TestCase):
    """The locus export sends the results query (buildApiUrl) to /export, every row rather than a page.

    Building it from buildApiUrl keeps every filter, and the way each one is spelled, the same in
    the results and in the export.
    """

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_the_export_url_is_built_from_the_results_query(self):
        export_url = js_function_body(self.template, "function buildLocusExportUrl(")
        self.assertIn("new URL(buildApiUrl())", export_url)
        self.assertIn("url.pathname.replace(/\\/loci$/, '/export')", export_url)
        for name in ("page", "page_size"):
            self.assertIn("url.searchParams.delete('%s');" % name, export_url)
        self.assertIn("url.searchParams.set('format', format);", export_url)

    def test_the_export_sets_no_filter_of_its_own(self):
        export = js_function_body(self.template, "async function performExport(")
        self.assertIn("? buildLocusExportUrl(format) : buildGroupedExportUrl(format);", export)
        self.assertNotIn("params.set(", export)

    def test_an_export_error_is_reported_from_the_download_frame(self):
        # The download runs in a hidden iframe, so a failed export's error body loads there, and
        # its load handler (which a successful attachment download does not fire) shows it
        export = js_function_body(self.template, "async function performExport(")
        self.assertIn("iframe.onload = () => {", export)
        self.assertIn("if (text.trim()) alert('Export failed: ' + text.trim().slice(0, 500));", export)
        self.assertIn("        iframe.src = url;\n    } finally {", export)

    def test_the_server_parses_the_export_query_like_the_results_query(self):
        server = read_template("results_server.py")
        for route in ('@app.route("/api/v1/loci")\ndef get_loci():', '@app.route("/api/v1/export")\ndef export_loci():'):
            body = server[server.index(route):]
            body = body[:body.index("\n@app.")]
            self.assertIn("validate_params()", body, route)
            self.assertIn("build_api_query(params)", body, route)


class VariationClusterHelpTextTests(unittest.TestCase):
    """The size-diff help popup must not contradict the checkbox help below it."""

    def test_missing_cluster_is_described_as_unknown_not_as_absent(self):
        template = read_template("index_page_template.html")
        self.assertIn("Loci with no variation cluster size are kept whatever value you enter, "
                      "because there is no size to compare against", template)


class SwimDetailModalTests(unittest.TestCase):
    """The swim plot's detail modal: one set of note and tag buttons, and the latest dot wins.

    The note and tag buttons are the ones rendered in its User Annotations section, so they only
    exist once the locus's annotations are on the page.
    """

    def setUp(self):
        self.template = read_template("swim_plot_template.html")
        self.show_dot_detail = js_function_body(self.template, "function showDotDetail(item) {")

    def test_the_footer_holds_only_the_close_button(self):
        self.assertIn('<div class="scrolling content" id="dot-detail-content"></div>\n'
                      '    <div class="actions">\n'
                      '        <div class="ui cancel button">Close</div>\n'
                      '    </div>', self.template)

    def test_the_rendered_section_has_the_note_and_tag_buttons(self):
        render = js_function_body(self.template, "function renderDotDetailContent(")
        for button_class in ("swim-inline-edit-note", "swim-inline-delete-note", "swim-inline-add-note",
                             "swim-inline-add-tag"):
            self.assertIn(button_class, render)

    def test_a_stale_response_cannot_render_into_the_modal(self):
        # The newer-dot guard returns before the render, so a late response for a dot the user
        # has already navigated away from cannot fill the modal now showing another one.
        self.assertLess(self.show_dot_detail.index("if (seq !== swimDetailSeq) return;"),
                        self.show_dot_detail.index("renderDotDetailContent(item, detail, stats, params)"))

    def test_show_dot_detail_itself_records_the_item_the_modal_acts_on(self):
        self.assertTrue(self.show_dot_detail.startswith("function showDotDetail(item) {\n"
                                                        "    _currentModalItem = item;\n"))
        self.assertNotIn("showDotDetail = ", self.template)

    def test_a_failed_sample_outlier_stats_request_is_logged(self):
        # fetchJsonOrNull (header_template.html) logs the failure; the modal then leaves out the
        # "Sample is outlier at loci" section
        self.assertIn("const statsPromise = fetchJsonOrNull(statsUrl, 'Load sample outlier stats');",
                      self.show_dot_detail)


class SearchErrorMessageTests(unittest.TestCase):
    """A failed search names what failed: reaching the server, the server, or drawing the results.

    describeFetchError turns any TypeError into "Could not reach the server", which is what fetch()
    rejects with when the request cannot be made. A TypeError thrown while rendering the response
    is a bug in the page, so performSearch shows its own message for it instead.
    """

    def setUp(self):
        self.search = js_function_body(read_template("index_page_template.html"), "async function performSearch(")

    def test_the_render_step_is_marked_once_the_response_is_in(self):
        self.assertIn("        rendering = true;\n"
                      "        if (grouped) renderGroupedResults(data); else renderResults(data);", self.search)
        self.assertLess(self.search.index("const data = await response.json();"),
                        self.search.index("rendering = true;"))

    def test_a_render_error_shows_its_own_message(self):
        self.assertIn("$('#error-text').text(rendering\n"
                      "            ? `The results could not be displayed: ${error.message} (details in the browser console)`\n"
                      "            : describeFetchError(error));", self.search)
        self.assertIn("console.error('Search error:', error);", self.search)


class SharedNoteEditorTests(unittest.TestCase):
    """Both pages edit notes in the one multi-line note editor in note_editor.html.

    A textarea keeps a note's line breaks, which a window.prompt text field cannot hold. Only the
    two pages that edit notes include the editor; the QC pages do not.
    """

    def setUp(self):
        self.editor_template = read_template("note_editor.html")

    def test_the_editor_is_a_textarea_in_its_own_partial(self):
        self.assertIn('<div id="edit-note-modal" class="ui small modal">', self.editor_template)
        self.assertIn('<textarea id="note-textarea" rows="4"', self.editor_template)
        for name in ("index_page_template.html", "swim_plot_template.html", "header_template.html"):
            self.assertNotIn('id="edit-note-modal"', read_template(name), name)

    def test_only_the_pages_that_edit_notes_include_it(self):
        including = [name for name in sorted(os.listdir(CODE_DIRECTORY)) if name.endswith(".html")
                     and '{% include "note_editor.html" %}' in read_template(name)]
        self.assertEqual(["index_page_template.html", "swim_plot_template.html"], including)

    def test_both_pages_open_it(self):
        self.assertIn("openNoteEditor(locusId, refreshAfterNoteChange);", read_template("index_page_template.html"))
        swim = read_template("swim_plot_template.html")
        self.assertIn("openNoteEditor($(this).data('locus-id'), refreshSwimModal);", swim)
        self.assertNotIn("prompt('Enter note:'", swim)

    def test_it_opens_over_the_swim_plot_detail_modal(self):
        self.assertIn("allowMultiple: true,", js_function_body(self.editor_template, "async function openNoteEditor("))

    def test_every_open_sets_the_editor_and_no_close_resets_it(self):
        # A reset when the hide animation ends could land after a quick reopen had already filled
        # the editor with the stored note, leaving it empty over that note.
        editor = js_function_body(self.editor_template, "async function openNoteEditor(")
        self.assertIn("$textarea.val(data.note.note_text);", editor)
        self.assertNotIn("onHidden: function", editor)

    def test_save_and_delete_tell_the_page_to_repaint(self):
        editor = js_function_body(self.editor_template, "async function openNoteEditor(")
        self.assertIn("noteEditorOnNoteChanged = onNoteChanged;", editor)
        self.assertIn("noteEditorOnNoteChanged(currentLocusId);", editor)
        delete_handler = self.editor_template[self.editor_template.index("$(document).on('click', '#delete-note-btn'"):]
        self.assertIn("noteEditorOnNoteChanged(locusId);", delete_handler[:delete_handler.index("\n});\n")])


class NoteEditorLoadFailureTests(unittest.TestCase):
    """The note editor must not open empty when the stored note could not be loaded.

    It fetches the note before opening the modal, so it has no click-while-loading window, but
    treating a non-ok response or a thrown request as "no note" would open an empty textarea, and
    saving that would replace the stored note with whatever was typed.
    """

    def setUp(self):
        self.handler = js_function_body(read_template("note_editor.html"), "async function openNoteEditor(")

    def test_a_failed_response_is_not_read_as_an_absent_note(self):
        self.assertIn("if (!response.ok) throw await errorFromResponse(response);", self.handler)

    def test_the_editor_stays_closed_when_the_existing_note_could_not_be_loaded(self):
        catch_block = self.handler[self.handler.index("} catch (e) {"):]
        catch_block = catch_block[:catch_block.index("\n    }")]
        self.assertIn("Failed to load the existing note: ", catch_block)
        self.assertIn("return;", catch_block)
        self.assertNotIn("$textarea.val('')", catch_block)

    def test_the_modal_is_only_shown_after_the_fetch_block(self):
        self.assertLess(self.handler.index("} catch (e) {"), self.handler.index(".modal('show')"))


class VariationClusterDetailPanelTests(unittest.TestCase):
    """The detail panel's variation-cluster line must show the same three states as the row.

    A locus whose variation cluster was computed and then filtered out carries a filter reason
    and a NULL size difference, so gating the whole line on the size difference would make the
    panel read as "no variation cluster" for exactly the loci whose VC cell shows a crossed-out icon.
    """

    def setUp(self):
        template = read_template("index_page_template.html")
        self.row = js_function_body(template, "function renderLocusRow(locus, rowIndex) {")
        self.detail = js_function_body(template, "function renderLocusDetail(detail, outlierType) {")
        self.helper = js_function_body(template, "function describeVariationCluster(locus) {")

    def test_one_helper_decides_the_state_for_both(self):
        for body in (self.row, self.detail):
            self.assertIn("const variationCluster = describeVariationCluster(locus);", body)
            self.assertNotIn("locus.VariationClusterSizeDiff", body)
            self.assertNotIn("locus.VariationClusterFilterReason", body)

    def test_the_helper_tells_the_three_states_apart(self):
        # A filtered-out cluster has a filter reason and a NULL size difference, so the reason is
        # looked at only when there is no size difference.
        self.assertIn("if (locus.VariationClusterSizeDiff != null && locus.VariationClusterSizeDiff !== '') {",
                      self.helper)
        self.assertIn("if (!Number.isFinite(sizeDiff)) return null;", self.helper)
        self.assertIn("sizeDiffText: (sizeDiff >= 0 ? '+' : '') + sizeDiff.toLocaleString('en-US') + 'bp'",
                      self.helper)
        self.assertIn("reasonTitle: VC_FILTER_REASON_TITLES[reason] || reason", self.helper)
        self.assertIn("reasonColor: VC_FILTER_REASON_COLORS[reason] || 'darkorange'", self.helper)

    def test_both_show_the_kept_and_the_filtered_out_state(self):
        for body in (self.row, self.detail):
            self.assertIn("if (variationCluster && variationCluster.sizeDiff != null) {", body)
            self.assertIn("} else if (variationCluster) {", body)
            self.assertIn('<i class="times circle outline icon" style="color: ${variationCluster.reasonColor};"', body)
        self.assertIn("Variation Cluster: ${variationCluster.sizeDiffText}", self.detail)


class AddTagEnterTests(unittest.TestCase):
    """Enter in the Edit Tags modal adds the highlighted menu item, not the typed prefix.

    Semantic UI commits the highlighted menu item in its own keydown handler on the dropdown
    element, so the page binds its Enter handler on the dropdown, after Semantic UI's, and reads the
    committed value (typing 'can' with 'candidate' highlighted adds 'candidate'). A key handler on
    the dropdown's search input would run before Semantic UI's, while nothing is committed yet.
    """

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_the_enter_handler_runs_after_semantic_uis(self):
        self.assertIn("$addTagDropdown.off('keydown.addtag').on('keydown.addtag', function(e) {", self.template)
        self.assertIsNone(re.search(r"input\.search'\)[^;]*\.on\('key", self.template))

    def test_the_typed_text_is_only_a_fallback(self):
        self.assertIn("const tag = (dropdownValue || searchValue || '').trim();", self.template)


class SharedRenderingHelperTests(unittest.TestCase):
    """Markup and rules that several places show are written once, in one helper each."""

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_both_expandable_texts_use_one_helper(self):
        # The .expand-phenotype / .collapse-phenotype handlers look for the classes it writes
        helper = js_function_body(self.template, "function renderExpandableText(text, maxLength) {")
        for name in ("phenotype-short", "expand-phenotype", "phenotype-full", "collapse-phenotype"):
            self.assertIn('class="%s"' % name, helper)
            self.assertEqual(1, self.template.count('class="%s"' % name), name)
        self.assertIn("renderExpandableText(locus.GeneTableLLMPhenotypeSummary, 50)",
                      js_function_body(self.template, "function renderLocusRow(locus, rowIndex) {"))
        self.assertIn("renderExpandableText(s.phenotype_description, 150)",
                      js_function_body(self.template, "function renderLocusDetail(detail, outlierType) {"))

    def test_both_summary_charts_use_one_builder(self):
        self.assertEqual(1, self.template.count("new Chart("))
        self.assertIn("new Chart(", js_function_body(self.template, "function renderSummaryBarChart("))
        for call in ("renderSummaryBarChart('motif-chart', 'motif-chart-container', labels, data, {",
                     "renderSummaryBarChart('gene-region-chart', 'gene-region-chart-container', labels, data, {"):
            self.assertIn(call, self.template)

    def test_the_swim_dots_and_legend_counts_use_one_status_rule(self):
        swim = read_template("swim_plot_template.html")
        self.assertIn("const colors = STATUS_COLORS[statusCategory(item.affected_status)];", swim)
        self.assertIn("data.forEach(item => { counts[statusCategory(item.affected_status)]++; });", swim)
        self.assertEqual(1, swim.count("isAffectedStatus("))


class DetailRowOuterCellTests(unittest.TestCase):
    """Detail-row lookups target the row's own cell, never the hundreds of cells nested inside it.

    The rendered detail holds a per-sample sub-table plus the threshold, phenotype-score and
    affected-unsolved tables, so `.find('td')` on the row would match every one of those cells:
    refreshDetailPanel would assign the whole detail string to each, and the collapse handlers
    would wrap and animate each, firing their slideUp callback once per element.
    """

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_no_detail_row_lookup_uses_a_descendant_selector(self):
        # Comment lines may name the descendant selector; code lines may not.
        for number, line in enumerate(self.template.split("\n"), start=1):
            if line.strip().startswith("//"):
                continue
            self.assertNotIn("find('td')", line, "line %d selects descendant cells" % number)
            self.assertNotIn("$td.find('.slide-wrapper')", line,
                             "line %d selects descendant wrappers" % number)

    def test_every_detail_cell_lookup_targets_a_direct_child(self):
        # Collapsing one row (collapseLocusDetailRow, shared by the expand toggle and the Collapse
        # details button), the loading row's error message, the detail render (loadLocusDetailInto,
        # shared by the expand handler and the post-annotation refresh) and collapse all, plus the
        # two colspans (a locus detail row's and a grouped row's) that count the cells of the row
        # they open under.
        self.assertEqual(6, self.template.count(".children('td')"))
        self.assertEqual(2, self.template.count("colspan=\"${$row.children('td').length}\""))
        # One slide wrapper per row: collapse one, collapse all.
        self.assertEqual(2, self.template.count(".children('.slide-wrapper')"))

    def test_both_ways_of_closing_one_panel_use_the_shared_collapse(self):
        self.assertIn("    if ($existingDetail.length) {\n        collapseLocusDetailRow($existingDetail);",
                      self.template)
        self.assertIn("$(document).on('click', '.collapse-detail-btn', function() {\n"
                      "    collapseLocusDetailRow($(this).closest('.locus-detail-row'));\n});", self.template)

    def test_the_wrapped_cell_is_the_one_that_is_animated(self):
        self.assertIn("$td.wrapInner('<div class=\"slide-wrapper\"></div>');\n"
                      "    $td.children('.slide-wrapper').slideUp(200, function() {", self.template)


class ResultsTableColumnCountTests(unittest.TestCase):
    """The results table's colspans follow its columns instead of a hard-coded count.

    The locus table's column count changes whenever a column is added, and a grouped view shows a
    header of its own, so every colspan that spans the results table is computed. The per-sample
    sub-table's colspans use its own ${colCount}.
    """

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_the_locus_header_and_the_locus_row_have_the_same_cells(self):
        header = self.template[self.template.index('<tr id="locus-header-row">'):]
        header = header[:header.index("</tr>")]
        row = js_function_body(self.template, "function renderLocusRow(locus, rowIndex) {")
        row = row[row.index('return `<tr class="locus-row"'):]
        row = row[:row.index("</tr>`;")]
        self.assertGreater(header.count("<th "), 0)
        # selectedAllelesCell, the grouped view's extra cell, is interpolated rather than literal
        self.assertEqual(header.count("<th "), row.count("<td "))

    def test_every_colspan_is_one_of_the_computed_counts(self):
        self.assertEqual({"${colCount}", "${resultsTableColumnCount()}", "${$row.children('td').length}"},
                         set(re.findall(r'colspan="([^"]+)"', self.template)))

    def test_the_count_follows_whichever_header_row_is_shown(self):
        helper = js_function_body(self.template, "function resultsTableColumnCount(")
        self.assertIn("$('#locus-header-row') : $('#group-header-row')", helper)
        self.assertIn("$row.children('th').length", helper)


class FamilyCountColumnTests(unittest.TestCase):
    """The #Fam column shows NumAffectedUnsolvedFamiliesAboveUnaffected and sorts by family_count."""

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_the_header_sorts_by_family_count_next_to_aff(self):
        self.assertIn('#Aff<span class="sort-icon">&#8597;</span></th>\n'
                      '                                <th class="col-fam-unsolved col-numeric" data-sort="family_count" '
                      'title="Distinct affected, unsolved families', self.template)

    def test_the_row_cell_reads_the_family_count(self):
        row = js_function_body(self.template, "function renderLocusRow(locus, rowIndex) {")
        self.assertIn("const nFamUnsolved = locus.NumAffectedUnsolvedFamiliesAboveUnaffected || 0;", row)
        self.assertIn('<td class="col-aff-unsolved col-numeric" style="${affUnsolvedStyle}">${nAffUnsolved}</td>\n'
                      '        <td class="col-fam-unsolved col-numeric" style="${famUnsolvedStyle}">${nFamUnsolved}</td>', row)

    def test_the_column_has_a_width(self):
        self.assertIn("#results-table .col-fam-unsolved { width: 55px; }", read_template("header_template.html"))

    def test_the_server_accepts_the_sort(self):
        self.assertIn('"family_count": ("NumAffectedUnsolvedFamiliesAboveUnaffected_{ot}", "DESC")',
                      read_template("results_server.py"))


class MinDeltaPopFilterTests(unittest.TestCase):
    """The Min Δpop input (min_repeats_above_population) goes through every state function.

    A control missing from any one of them is silently lost: on reload (readStateFromUrl only keeps
    DEFAULT_STATE keys), on Search (readStateFromUI), in the request and the export (buildApiUrl),
    on Back or a chip removal (applyStateToUI), or from the active filter chips.
    """

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_the_input_exists(self):
        self.assertIn('<input type="number" name="min_repeats_above_population"', self.template)

    def test_every_state_function_handles_it(self):
        self.assertIn("    min_repeats_above_population: '',",
                      js_object_body(self.template, "const DEFAULT_STATE = {"))
        # applyStateToUI, readStateFromUI and buildApiUrl take it from the shared list
        # (FilterControlListTests)
        self.assertIn("min_repeats_above_population", js_string_list(self.template, "TEXT_FILTER_INPUTS"))
        self.assertIn("min_repeats_above_population: 'Min Δpop'",
                      js_function_body(self.template, "function updateActiveFiltersBar("))

    def test_it_makes_the_population_metric_count(self):
        # The filter compares against the 99th % / max choice, so 'max' has to reach the API with it.
        self.assertIn("currentState.min_repeats_above_population",
                      js_function_body(self.template, "function populationMetricAffectsResults("))
        self.assertIn("populationMetricAffectsResults() && currentState.population_metric === 'max'",
                      js_function_body(self.template, "function buildApiUrl("))

    def test_the_delta_pop_column_compares_against_the_cohorts_the_server_filters_on(self):
        # The filter is a minimum on the Δpop column, so the column has to use the population
        # metric and the cohorts of results_server.population_datasets_compared.
        row = js_function_body(self.template, "function renderLocusRow(locus, rowIndex) {")
        self.assertIn("const popMetric = currentState.population_metric === 'max' ? 'MaxAllele' : '99thPercentile';", row)
        self.assertIn("currentState.require_above_population === 'short-read' ? ['TenK10K']", row)
        self.assertIn("currentState.require_above_population === 'long-read' ? ['HPRC256', 'AoU1027']", row)
        self.assertIn(": POPULATION_COHORTS;", row)
        self.assertIn("locus[`${p}_${popMetric}`]", row)
        # The full cohort list is pinned to the server's in PopulationCohortTests
        self.assertIn('"short-read": ("TenK10K",),', read_template("results_server.py"))

    def test_a_url_with_a_parameter_the_page_does_not_use_is_normalized_in_place(self):
        # readStateFromUrl keeps only DEFAULT_STATE keys, and then always rewrites the arrival
        # hash in place, so a URL carrying a parameter the page does not use (e.g. min_pli, which
        # only API callers send) neither keeps it nor traps Back.
        read_state = js_function_body(self.template, "function readStateFromUrl(")
        self.assertIn("if (key in DEFAULT_STATE) {", read_state)
        self.assertIn("    replaceStateInUrl();\n\n    applyStateToUI();", read_state)
        self.assertNotIn("min_pli", js_object_body(self.template, "const DEFAULT_STATE = {"))
        self.assertIn('request.args.get("min_pli")', read_template("results_server.py"))

    def test_the_metric_stays_in_the_url_because_the_column_uses_it(self):
        state_hash = js_function_body(self.template, "function buildStateHash(")
        self.assertNotIn("key === 'population_metric'", state_hash)


class FilterControlListTests(unittest.TestCase):
    """applyStateToUI, readStateFromUI and buildApiUrl take the filter controls from one shared list
    per kind.

    A control left out of a list is silently lost on reload, on Search or from the request, so each
    list has to name exactly the page's inputs of that kind.
    """

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_the_lists_name_exactly_the_pages_text_and_checkbox_filters(self):
        # The population filter field's inputs come from filter_control_macros.html
        page = self.template + str(filter_control_macros().population_filter_field("", "name"))
        inputs = re.findall(r'<input type="(text|number|checkbox)" name="(\w+)"', page)
        self.assertEqual(sorted(name for kind, name in inputs if kind in ("text", "number")),
                         sorted(js_string_list(self.template, "TEXT_FILTER_INPUTS")))
        # The sort checkboxes belong to the sort pills, and the gene-region ones to GENE_REGION_FILTERS
        self.assertEqual(sorted(name for kind, name in inputs if kind == "checkbox" and name != "sort_by"
                                and not re.match(r"(exclude_)?gene_region_", name)),
                         sorted(js_string_list(self.template, "CHECKBOX_FILTERS")))

    def test_each_name_is_a_state_key(self):
        default_state = js_object_body(self.template, "const DEFAULT_STATE = {")
        for name in js_string_list(self.template, "TEXT_FILTER_INPUTS"):
            self.assertIn("    %s: ''," % name, default_state)
        for name in js_string_list(self.template, "CHECKBOX_FILTERS"):
            self.assertIn("    %s: false," % name, default_state)

    def test_both_state_functions_loop_over_the_lists(self):
        for signature in ("function applyStateToUI(", "function readStateFromUI("):
            body = js_function_body(self.template, signature)
            for loop in ("for (const name of TEXT_FILTER_INPUTS) {", "for (const name of CHECKBOX_FILTERS) {",
                         "for (const [key, namePrefix] of GENE_REGION_FILTERS) {"):
                self.assertIn(loop, body, signature)

    def test_the_request_sends_every_filter_in_the_lists(self):
        build = js_function_body(self.template, "function buildApiUrl(")
        self.assertIn("...TEXT_FILTER_INPUTS, 'sample_id_keyword', ...GENE_REGION_FILTERS.map(([key]) => key), 'sort_by', 'tag',",
                      build)
        self.assertIn("for (const name of stringParams) {\n        if (currentState[name]) params.set(name, currentState[name]);",
                      build)
        self.assertIn("for (const name of CHECKBOX_FILTERS) {", build)


class GroupResultsByTests(unittest.TestCase):
    """The Group Results By controls go through every state function and reach the grouped API."""

    def setUp(self):
        self.template = read_template("index_page_template.html")
        self.server = read_template("results_server.py")

    def test_every_grouped_column_has_a_width(self):
        # In the table's fixed layout an unsized column only gets the space the sized ones leave,
        # which leaves Phenotype zero pixels wide at 1280px
        header = read_template("header_template.html")
        for column_class in set(re.findall(r"'(col-group-[a-z-]+)", self.template)):
            self.assertRegex(header, r"#results-table \." + column_class + r"[ ,][^}]*width:", column_class)

    def test_removing_the_group_sort_chip_also_restores_the_default_direction(self):
        # Otherwise a leftover ASC from 'Sample / family id' sorts the default # of loci fewest first
        start = self.template.index("$(document).on('click', '.filter-chip .delete', function() {")
        handler = self.template[start:self.template.index("});", start)]
        self.assertIn("if (key === 'group_sort') {", handler)
        self.assertIn("currentState.group_sort_direction = GROUP_SORT_NATURAL_DIRECTIONS[DEFAULT_STATE.group_sort] || 'desc';",
                      handler)

    def test_the_nested_locus_header_is_not_sticky(self):
        # The grouped row's locus table sits inside #results-table, whose header cells are sticky
        self.assertIn("#results-table .group-loci-table thead th { position: static; }",
                      read_template("header_template.html"))

    def test_the_affected_quick_filter_uses_the_samples_own_status(self):
        # Rows are colored by the #Aff rule, which counts unknown-status samples as affected-unsolved;
        # the Affected button must still show only samples whose status is affected
        self.assertIn("if (isAffectedStatus(s.affected_status)) rowClass += ' sample-status-affected';", self.template)
        self.assertIn("if (filter === 'affected') $row.toggle($row.hasClass('sample-status-affected'));", self.template)

    def test_the_controls_exist(self):
        for value in ("locus", "sample", "family"):
            self.assertIn('<input type="radio" name="group_by" value="%s"' % value, self.template)
        for value in ("asc", "desc"):
            self.assertIn('<input type="radio" name="group_sort_direction" value="%s"' % value, self.template)
        self.assertIn('<div id="group-sort-dropdown" class="ui compact selection dropdown">', self.template)

    def test_the_sort_menu_offers_exactly_the_server_sorts(self):
        menu = self.template[self.template.index('<div id="group-sort-dropdown"'):]
        menu = menu[:menu.index("</span>")]
        server_sorts = self.server[self.server.index("GROUPED_SORTS = {"):]
        server_sorts = server_sorts[:server_sorts.index("\n}\n")]
        self.assertEqual(set(re.findall(r'^    "(\w+)":', server_sorts, re.MULTILINE)),
                         set(re.findall(r'class="item" data-value="(\w+)"', menu)))
        # The id sorts A to Z by default, as the page's GROUP_SORT_NATURAL_DIRECTIONS assumes.
        self.assertIn('"group_id": "asc",', server_sorts)
        self.assertIn("const GROUP_SORT_NATURAL_DIRECTIONS = {group_id: 'asc'};", self.template)

    def test_every_state_function_handles_them(self):
        default_state = js_object_body(self.template, "const DEFAULT_STATE = {")
        for line in ("group_by: 'locus',", "group_sort: 'num_loci',", "group_sort_direction: 'desc',"):
            self.assertIn(line, default_state)
        apply_state = js_function_body(self.template, "function applyStateToUI(")
        self.assertIn("'require_above_population', 'group_by']", apply_state)
        self.assertIn("setGroupSortDropdown(currentState.group_sort);", apply_state)
        self.assertIn("inputsWithValue('group_sort_direction', currentState.group_sort_direction)", apply_state)
        read_state = js_function_body(self.template, "function readStateFromUI(")
        for name in ("group_by", "group_sort", "group_sort_direction"):
            self.assertIn("currentState.%s = " % name, read_state)
        labels = js_function_body(self.template, "function updateActiveFiltersBar(")
        self.assertIn("group_by: 'Group by', group_sort: 'Group sort', group_sort_direction: 'Group sort direction',",
                      labels)

    def test_the_sort_is_left_out_of_the_url_and_the_chips_in_the_locus_view(self):
        skip = "(key === 'group_sort' || key === 'group_sort_direction') && currentState.group_by === 'locus'"
        self.assertIn(skip, js_function_body(self.template, "function buildStateHash("))
        self.assertIn(skip, js_function_body(self.template, "function updateActiveFiltersBar("))

    def test_search_and_export_use_the_grouped_routes(self):
        search = js_function_body(self.template, "async function performSearch(")
        self.assertIn("fetchWithTimeout(grouped ? buildGroupedApiUrl() : buildApiUrl())", search)
        self.assertIn("if (grouped) renderGroupedResults(data); else renderResults(data);", search)
        export = js_function_body(self.template, "async function performExport(")
        # The export menu holds the formats of the table on screen, so its grouping wins over
        # an unsearched change of the radios.
        self.assertIn("currentState.group_by = savedState.group_by;", export)
        self.assertIn(": buildGroupedExportUrl(format);", export)
        self.assertIn("'/outlier_samples'", js_function_body(self.template, "function buildGroupedApiUrl("))
        self.assertIn("'/outlier_samples/export'", js_function_body(self.template, "function buildGroupedExportUrl("))
        self.assertIn('@app.route("/api/v1/outlier_samples")', self.server)
        self.assertIn('@app.route("/api/v1/outlier_samples/export")', self.server)
        self.assertIn('GROUPED_EXPORT_FORMATS = ("groups_tsv", "group_loci_tsv")', self.server)
        for export_format in ("groups_tsv", "group_loci_tsv"):
            self.assertIn('data-value="%s"' % export_format, self.template)

    def test_a_grouped_row_expands_to_its_own_loci(self):
        loci_url = js_function_body(self.template, "function buildGroupLociApiUrl(")
        for name in ("affected_unsolved_outlier_family_id", "affected_unsolved_outlier_sample_id"):
            self.assertIn("'%s'" % name, loci_url)
            self.assertIn('"%s"' % name, self.server)
        self.assertIn("for (const name of SAMPLE_LEVEL_FILTERS) url.searchParams.delete(name);", loci_url)
        self.assertIn("const SAMPLE_LEVEL_FILTERS = ['sample_id_keyword', 'sample_id_like', 'phenotype_keyword'];",
                      self.template)
        self.assertIn('SAMPLE_LEVEL_FILTERS = ("sample_id_keyword", "sample_id_like", "phenotype_keyword")',
                      self.server)

    def test_the_expand_handler_is_bound_to_table_cells_only(self):
        # The grouped header row has a th with the same class, which has no row to expand.
        self.assertIn("$(document).on('click', 'td.col-group-expand'", self.template)
        self.assertNotIn("$(document).on('click', '.col-group-expand'", self.template)

    def test_the_sorted_group_column_uses_the_one_styled_indicator_class(self):
        self.assertIn("`${cls} sorted-active`",
                      js_function_body(self.template, "function renderGroupedResults(data) {"))
        self.assertIn("removeClass('sorted-active')",
                      js_function_body(self.template, "function renderGroupLociTable(groupId, data) {"))


class PopulationGateOptionTests(unittest.TestCase):
    """The "All" population filter and the server's --initial-require-above-population start value."""

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_both_pages_offer_the_all_cohorts_filter_the_server_accepts(self):
        # Both pages render the field from filter_control_macros.html with their own name prefix; the
        # swim plot names its checkbox by id, as it reads its checkboxes
        field = filter_control_macros().population_filter_field
        for name, prefix, attribute in (("index_page_template.html", "", "name"), ("swim_plot_template.html", "swim_", "id")):
            self.assertIn("{{ population_filter_field('%s', '%s') | indent(" % (prefix, attribute), read_template(name), name)
            rendered = str(field(prefix, attribute))
            self.assertIn('<input type="radio" name="%srequire_above_population" value="all">' % prefix, rendered)
            self.assertIn('<input type="checkbox" %s="%sinclude_loci_without_population_data">' % (attribute, prefix), rendered)
        self.assertIn('"all": POPULATION_DATASETS,', read_template("results_server.py"))

    def test_the_server_hands_the_page_its_start_value(self):
        self.assertIn("const INITIAL_REQUIRE_ABOVE_POPULATION = {{ initial_require_above_population | tojson }};",
                      self.template)
        self.assertIn('initial_require_above_population=app.config.get("INITIAL_REQUIRE_ABOVE_POPULATION", ""),',
                      read_template("results_server.py"))

    def test_it_only_applies_to_a_url_with_no_fragment_at_all(self):
        startup = self.template[self.template.index("Promise.all([tagDropdownReady, sampleIdDropdownReady])"):]
        startup = startup[:startup.index("performSearch();")]
        self.assertLess(startup.index("readStateFromUrl();"), startup.index("INITIAL_REQUIRE_ABOVE_POPULATION"))
        self.assertIn("if (!window.location.href.includes('#') && INITIAL_REQUIRE_ABOVE_POPULATION) {", startup)

    def test_the_default_state_keeps_a_bare_fragment_in_the_url(self):
        # Otherwise removing the start filter leaves a URL with no '#', and a reload sets it again.
        for signature in ("function writeStateToUrl(", "function replaceStateInUrl("):
            body = js_function_body(self.template, signature)
            self.assertIn("hash || '#'", body, signature)
            self.assertNotIn("window.location.pathname", body, signature)


class TRExplorerLinkTests(unittest.TestCase):
    """Every TRExplorer locus link is built by header_template.html's buildTRExplorerUrl."""

    def test_the_header_builds_the_url_with_the_chosen_columns(self):
        helper = js_function_body(read_template("header_template.html"), "function buildTRExplorerUrl(locus) {")
        self.assertIn("&showColumns=${encodeURIComponent(TREXPLORER_LINK_COLUMNS.join(','))}", helper)

    def test_both_pages_use_it_and_build_no_url_of_their_own(self):
        for name in ("index_page_template.html", "swim_plot_template.html"):
            template = read_template(name)
            self.assertIn("buildTRExplorerUrl(locus)", template, name)
            self.assertNotIn("trexplorer.broadinstitute.org/#", template, name)
            self.assertIn('{% include "header_template.html" %}', template, name)


class SwimPlotMarkerTests(unittest.TestCase):
    """The swim plot marks true/candidate tags with a filled triangle, and the legend shows the same marker."""

    def test_the_point_style_and_the_legend_agree(self):
        template = read_template("swim_plot_template.html")
        self.assertIn("style = 'triangle';", template)
        self.assertIn('<span style="margin-left: 8px;">&#9650;</span> True/Candidate', template)

    def test_only_a_whole_true_or_candidate_tag_gets_the_triangle(self):
        # item.tags is the locus's tag list, so includes() compares whole tags, and a tag such as
        # 'not a candidate' keeps the circle.
        template = read_template("swim_plot_template.html")
        self.assertIn("const tags = (item.tags || []).map(tag => String(tag).trim().toLowerCase());", template)
        self.assertIn("if (tags.includes('true') || tags.includes('candidate')) {", template)

    def test_the_cap_note_matches_the_server_limit(self):
        self.assertIn("const SWIM_PLOT_MAX_ALLELES_PER_MOTIF_SIZE = 500;", read_template("swim_plot_template.html"))
        self.assertIn("ORDER BY allele_size DESC, LocusId, sample_id LIMIT 500", read_template("results_server.py"))


class SharedMenuLoaderTests(unittest.TestCase):
    """Both pages fill their tag and sample-ID menus with header_template.html's loaders.

    The loaders go through fetchJsonOrNull, which treats an error status as a failed request, so
    a failed request leaves a menu and the cached tag list as they were instead of emptying them.
    """

    def setUp(self):
        self.header = read_template("header_template.html")

    def test_an_error_status_is_a_failed_request(self):
        loader = js_function_body(self.header, "async function fetchJsonOrNull(")
        self.assertIn("if (!response.ok) throw await errorFromResponse(response);", loader)
        self.assertIn("return null;", loader)

    def test_one_tag_request_fills_the_menu_and_the_cached_list(self):
        refresh = js_function_body(self.header, "async function refreshTagFilterDropdown(")
        self.assertIn("const data = await fetchJsonOrNull(`${ANNOTATIONS_API_URL}/tags`, 'Load tags');\n"
                      "    if (!data) return;", refresh)
        self.assertIn("window._allTagsList = data.tags || [];", refresh)
        # The filter menu lists the same tags the add-tag hints offer
        self.assertIn("setupDropdownMenu($dropdown, window._allTagsList.map(", refresh)

    def test_both_pages_use_the_shared_loaders(self):
        for name, prefix in (("index_page_template.html", ""), ("swim_plot_template.html", "swim-")):
            template = read_template(name)
            self.assertIn("refreshTagFilterDropdown($('#%stag-filter'))" % prefix, template, name)
            self.assertIn("fillSampleIdDropdown($('#%ssample-id-filter'), schema);" % prefix, template, name)
            self.assertNotIn("/annotations/tags", template, name)

    def test_the_results_page_requests_the_schema_once(self):
        # The one response fills the sample-ID menu and the loci / samples count
        template = read_template("index_page_template.html")
        self.assertEqual(1, template.count("/schema`"))
        self.assertIn("$('#total-loci-label').text(counts);", template)


class SharedPageHelperTests(unittest.TestCase):
    """The results page and the swim plot share the UCSC link, the outlier-type map and the
    gene-region checkbox wiring through header_template.html, so the two copies cannot drift."""

    PAGES = ("index_page_template.html", "swim_plot_template.html")

    def test_the_helpers_are_defined_in_the_header(self):
        header = read_template("header_template.html")
        self.assertIn("function buildUcscUrl(locus) {", header)
        self.assertIn("const OUTLIER_TYPE_MAP = {all: 'AllAlleles', short: 'ShortAlleles', hemi: 'HemizygousAlleles'};",
                      header)
        self.assertIn("function linkGeneRegionIncludeAndExcludeCheckboxes(namePrefix) {", header)
        # The JavaScript map follows the server's
        self.assertIn('OUTLIER_TYPE_MAP = {"all": "AllAlleles", "short": "ShortAlleles", "hemi": "HemizygousAlleles"}',
                      read_template("results_server.py"))

    def test_both_pages_use_them(self):
        for name, prefix in zip(self.PAGES, ("''", "'swim_'")):
            template = read_template(name)
            self.assertIn("const ucscUrl = buildUcscUrl(locus);", template, name)
            self.assertIn("OUTLIER_TYPE_MAP[", template, name)
            self.assertIn("linkGeneRegionIncludeAndExcludeCheckboxes(%s);" % prefix, template, name)
            self.assertNotIn("genome.ucsc.edu", template, name)


class DotJitterTests(unittest.TestCase):
    """Chart dots are spread by a fixed offset per dot, so a redraw leaves every dot in place."""

    def test_the_swim_plot_offset_comes_from_the_locus_and_sample(self):
        self.assertIn("y: categoryToY(item.motif_category) + jitterFromKey(`${item.LocusId}\\t${item.sample_id}`, 0.6),",
                      read_template("swim_plot_template.html"))

    def test_the_qc_box_plot_offset_comes_from_the_point(self):
        self.assertIn("x: jitterFromKey(JSON.stringify(point), 0.6)",
                      js_function_body(read_template("qc_shared_js.html"), "function createSingleBinBoxPlotChart("))


class AnnotationRefreshTests(unittest.TestCase):
    """A note or tag edit requests the locus's annotations once, and a tag edit leaves the detail
    panel alone, since the panel shows the note but not the tags."""

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_a_note_change_repaints_the_badges_once_and_the_open_panel(self):
        self.assertEqual("function refreshAfterNoteChange(locusId) {\n"
                         "    refreshRowAnnotations(locusId);\n"
                         "    refreshDetailPanel(locusId);\n"
                         "}\n",
                         js_function_body(self.template, "function refreshAfterNoteChange("))
        self.assertNotIn("refreshRowAnnotations", js_function_body(self.template, "async function refreshDetailPanel("))

    def test_the_detail_panel_shows_the_note_but_not_the_tags(self):
        detail = js_function_body(self.template, "function renderLocusDetail(detail, outlierType) {")
        self.assertIn("detail.user_annotations?.note", detail)
        self.assertNotIn("tags", detail)


class PopulationCohortTests(unittest.TestCase):
    """Both pages take the population cohorts from one list, header_template.html's POPULATION_COHORTS."""

    def test_the_list_matches_the_server_cohorts(self):
        self.assertIn("const POPULATION_COHORTS = ['HPRC256', 'AoU1027', 'TenK10K'];",
                      read_template("header_template.html"))
        server = read_template("results_server.py")
        self.assertIn('LONG_READ_POPULATION_DATASETS = ("HPRC256", "AoU1027")', server)
        self.assertIn('POPULATION_DATASETS = LONG_READ_POPULATION_DATASETS + ("TenK10K",)', server)

    def test_both_detail_views_iterate_over_it(self):
        index = read_template("index_page_template.html")
        self.assertIn("const popEntries = POPULATION_COHORTS.map(name => ({",
                      js_function_body(index, "function renderLocusDetail(detail, outlierType) {"))
        swim = read_template("swim_plot_template.html")
        self.assertIn("POPULATION_COHORTS.filter(k => popThresh[k + '_MaxAllele'] != null).forEach(k => {",
                      js_function_body(swim, "function renderDotDetailContent("))


class AnnotationBadgeTests(unittest.TestCase):
    """A row repainted after a note or tag edit gets the same badges as a freshly loaded row."""

    def test_both_paths_use_the_one_badge_builder(self):
        template = read_template("index_page_template.html")
        self.assertEqual(1, template.count('<span class="badge badge-tag">'))
        self.assertIn("renderAnnotationBadges(locus.tags || [], locus.has_note)",
                      js_function_body(template, "function renderLocusRow(locus, rowIndex) {"))
        self.assertIn("renderAnnotationBadges(ann.tags || [], ann.note != null)",
                      js_function_body(template, "async function refreshRowAnnotations("))

    def test_every_row_listing_the_locus_is_repainted_on_its_own(self):
        # In a grouped view the same locus can be open under several groups, and only some of their
        # rows may have a badges row yet, so each row decides on its own whether to update or create one
        refresh = js_function_body(read_template("index_page_template.html"), "async function refreshRowAnnotations(")
        each_row = refresh[refresh.index("$(`tr.locus-row[data-locus-id=\"${locusId}\"]`).each(function() {"):]
        self.assertLess(each_row.index("const $badgesRow = $locusCell.find('.badges-row');"),
                        each_row.index("if ($badgesRow.length) {"))
        self.assertIn("$locusCell.find('.cell-main').after(`<div class=\"cell-sub badges-row\">${newContent}</div>`);",
                      each_row)


class CdnIntegrityTests(unittest.TestCase):
    """Every CDN file the pages load names an exact version and carries a Subresource Integrity hash."""

    CDN_TAG_PATTERN = r"<(?:script|link)\b[^>]*\b(?:src|href)=\"https?://[^>]*>"

    def test_every_cdn_tag_in_the_header_is_pinned_and_hashed(self):
        tags = re.findall(self.CDN_TAG_PATTERN, read_template("header_template.html"))
        self.assertEqual(4, len(tags))
        for tag in tags:
            self.assertRegex(tag, r'integrity="sha384-[A-Za-z0-9+/]{64}"')
            self.assertIn('crossorigin="anonymous"', tag)
            self.assertRegex(tag, r"[@/-]\d+\.\d+\.\d+[/.]")

    def test_no_other_template_loads_a_cdn_file(self):
        for name in sorted(os.listdir(CODE_DIRECTORY)):
            if name.endswith(".html") and name != "header_template.html":
                self.assertNotRegex(read_template(name), self.CDN_TAG_PATTERN, name)


class PagePastTheLastPageTests(unittest.TestCase):
    """A search that lands past the last page (its only row left the results, or a URL names a page
    that no longer exists) shows the last page instead of an empty one.

    The clamped page replaces the URL's history entry rather than adding one, so Back does not
    return to the out-of-range page, which would clamp and push again.
    """

    def setUp(self):
        self.search = js_function_body(read_template("index_page_template.html"), "async function performSearch(")

    def test_the_search_retries_on_the_last_page_before_rendering(self):
        self.assertIn("if (data.total_pages > 0 && data.page > data.total_pages) {\n"
                      "            currentState.page = data.total_pages;\n", self.search)
        self.assertLess(self.search.index("data.page > data.total_pages"), self.search.index("renderResults(data)"))
        # The server reports the page count the check compares against
        self.assertIn('"total_pages": total_pages,', read_template("results_server.py"))

    def test_the_clamped_page_replaces_the_history_entry_before_the_retry(self):
        clamp = self.search[self.search.index("currentState.page = data.total_pages;"):]
        self.assertLess(clamp.index("replaceStateInUrl();"), clamp.index("return performSearch();"))


class TypedInDropdownValueEscapingTests(unittest.TestCase):
    """Every dropdown that accepts typed-in values escapes them wherever Semantic UI writes them as HTML.

    The results page restores #sample-id-filter from the URL hash, so without the escaping a link
    such as #sample_id_keyword=<img src=x onerror=...> would run its script on the page. Semantic UI
    writes a typed-in value into the "Add ..." menu item (templates.addition) and into the selection
    label (templates.label), and with preserveHTML: false it hands the label template plain text.
    """

    def test_the_shared_settings_escape_the_addition_and_the_label(self):
        settings = js_object_body(read_template("qc_shared_js.html"), "const ESCAPING_ADDITIONS_DROPDOWN_SETTINGS = {")
        self.assertIn("preserveHTML: false,", settings)
        self.assertIn("addition: choice => escapeHtml(choice),", settings)
        self.assertIn("label: (value, text) => escapeHtml(text) + '<i class=\"delete icon\"></i>',", settings)

    def test_every_dropdown_that_accepts_typed_in_values_uses_them(self):
        dropdowns_with_additions = []
        for name in sorted(os.listdir(CODE_DIRECTORY)):
            if not name.endswith(".html"):
                continue
            for match in re.finditer(r"(\$\('#[\w-]+'\)|\$dropdown)\.dropdown\(\{(.*?)\n    \}\);", read_template(name), re.DOTALL):
                if "allowAdditions: true" in match.group(2):
                    dropdowns_with_additions.append(match.group(1))
                    self.assertIn("...ESCAPING_ADDITIONS_DROPDOWN_SETTINGS,", match.group(2), match.group(1))
        # header_template.html's setUpSampleIdFilterDropdown ($dropdown), and the add-tag dropdown
        self.assertEqual(["$dropdown", "$('#add-tag-dropdown')"], dropdowns_with_additions)
        self.assertIn("function setUpSampleIdFilterDropdown($dropdown) {\n    $dropdown.dropdown({",
                      read_template("header_template.html"))
        self.assertIn("setUpSampleIdFilterDropdown($('#sample-id-filter'));", read_template("index_page_template.html"))
        self.assertIn("setUpSampleIdFilterDropdown($('#swim-sample-id-filter'));", read_template("swim_plot_template.html"))

    def test_the_url_restored_dropdown_is_one_of_them(self):
        apply_state = js_function_body(read_template("index_page_template.html"), "function applyStateToUI(")
        self.assertIn("$('#sample-id-filter').dropdown('set exactly', ids);", apply_state)

    def test_url_values_are_compared_rather_than_written_into_selectors(self):
        # A URL value such as `"]` written into a selector makes it invalid, and applyStateToUI would
        # throw before the page's first search
        template = read_template("index_page_template.html")
        helper = js_function_body(template, "function inputsWithValue(")
        self.assertIn("$(`input[name=\"${name}\"]`).filter(function() { return this.value === value; });", helper)
        apply_state = js_function_body(template, "function applyStateToUI(")
        self.assertIn("inputsWithValue(name, currentState[name] || '').parent().checkbox('check');", apply_state)
        self.assertIn("inputsWithValue('population_metric', currentState.population_metric || '99th')", apply_state)

    def test_the_url_sort_list_is_compared_rather_than_written_into_a_selector(self):
        # performSearch marks the sorted headers before its try block, so a selector made invalid by
        # a URL sort field such as `size"]` would stop the search
        indicator = js_function_body(read_template("index_page_template.html"), "function syncSortHeaderIndicator(")
        self.assertIn("$('#results-table thead th[data-sort]')\n"
                      "        .filter(function() { return fields.includes(this.getAttribute('data-sort')); })\n"
                      "        .addClass('sorted-active');", indicator)


class ErrorResponseTests(unittest.TestCase):
    """Every failed API response is turned into an Error by header_template.html's errorFromResponse.

    It carries the server's "detail" or "error" text, and a response that is not JSON (e.g. a
    proxy's error page) gives "HTTP <status>".
    """

    def test_the_shared_helper_reads_the_body_with_a_fallback(self):
        helper = js_function_body(read_template("header_template.html"), "async function errorFromResponse(")
        self.assertIn("const body = await response.json().catch(() => ({}));", helper)
        self.assertIn("return new Error(body.detail || body.error || `HTTP ${response.status}`);", helper)

    def test_every_page_request_uses_it(self):
        templates = {name: read_template(name) for name in sorted(os.listdir(CODE_DIRECTORY)) if name.endswith(".html")}
        # The helper's own read is the only place a template reads an error body
        self.assertEqual(1, sum(template.count(".json().catch(") for template in templates.values()))
        for name, function_signature in (("index_page_template.html", "async function performSearch("),
                                         ("index_page_template.html", "async function loadLocusDetailInto("),
                                         ("swim_plot_template.html", "async function loadSwimPlot("),
                                         ("qc_shared_js.html", "function drawQCPageOnLoad(")):
            self.assertIn("if (!response.ok) throw await errorFromResponse(response);",
                          js_function_body(templates[name], function_signature), function_signature)


class DropdownValueReadTests(unittest.TestCase):
    """Every dropdown value is read through getDropdownValue, which undoes Semantic UI storing a
    double quote as '&quot;', so a tag a"b is saved and filtered as a"b."""

    def test_the_helper_restores_double_quotes(self):
        helper = js_function_body(read_template("qc_shared_js.html"), "function getDropdownValue(")
        self.assertIn("return String($dropdown.dropdown('get value') || '').replace(/&quot;/g, '\"');", helper)

    def test_it_is_the_only_code_that_reads_a_dropdown_value(self):
        templates = {name: read_template(name) for name in sorted(os.listdir(CODE_DIRECTORY)) if name.endswith(".html")}
        self.assertEqual(1, sum(template.count(".dropdown('get value')") for template in templates.values()))
        for name, reader in (("index_page_template.html", "getDropdownValue($('#tag-filter'))"),
                             ("index_page_template.html", "getDropdownValue($('#sample-id-filter'))"),
                             ("index_page_template.html", "getDropdownValue($addTagDropdown)"),
                             ("swim_plot_template.html", "getDropdownValue($('#swim-tag-filter'))"),
                             ("swim_plot_template.html", "getDropdownValue($('#swim-sample-id-filter'))"),
                             ("header_template.html", "const selected = getDropdownValue($dropdown);")):
            self.assertIn(reader, templates[name], name)


class SwimPlotMendelianGeneFilterTests(unittest.TestCase):
    """The swim plot sends the Mendelian-gene filter under the name the API parses, mendelian_genes_only."""

    def test_the_checkbox_sends_mendelian_genes_only(self):
        self.assertIn("if ($('#swim_mendelian_genes_only').is(':checked')) params.mendelian_genes_only = '1';",
                      read_template("swim_plot_template.html"))
        self.assertIn('"known_motifs_only", "mendelian_genes_only",', read_template("results_server.py"))


class NestedTableHeaderTests(unittest.TestCase):
    """Every table nested in #results-table overrides the main header's sticky positioning."""

    def test_every_nested_table_header_is_static(self):
        header = read_template("header_template.html")
        for table_class in ("group-loci-table", "sample-subtable", "phenotype-scores-table"):
            self.assertRegex(header, r"#results-table \.%s thead th \{[^}]*position: static;" % table_class)
        self.assertIn("#results-table .phenotype-scores-table thead th { position: static; background: transparent; }",
                      header)
        self.assertIn('<table class="phenotype-scores-table"', read_template("index_page_template.html"))


class LocusDetailLoadTests(unittest.TestCase):
    """The expand handler and the post-edit refresh load a locus's detail through one function, and
    every page builds the detail URL with header_template.html's buildLocusDetailUrl."""

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_both_paths_use_the_shared_loader(self):
        self.assertIn("await loadLocusDetailInto(locusId, $row);", self.template)
        self.assertIn("await loadLocusDetailInto(locusId, $rows);",
                      js_function_body(self.template, "async function refreshDetailPanel("))
        loader = js_function_body(self.template, "async function loadLocusDetailInto(")
        self.assertIn("fetchWithTimeout(buildLocusDetailUrl(locusId))", loader)
        self.assertIn("highlightGroupSamplesInDetail($(this), $detailRow);", loader)

    def test_every_detail_request_uses_the_shared_url(self):
        self.assertIn("function buildLocusDetailUrl(locusId) {", read_template("header_template.html"))
        # The shared loader and the known-locus popup; the swim plot's detail modal
        self.assertEqual(2, self.template.count("fetchWithTimeout(buildLocusDetailUrl(locusId))")
                         + self.template.count("fetch(buildLocusDetailUrl(locusId))"))
        self.assertIn("const detailUrl = buildLocusDetailUrl(item.LocusId);", read_template("swim_plot_template.html"))


class ApiBaseUrlTests(unittest.TestCase):
    """Every page calls the API through header_template.html's one API_BASE_URL."""

    def test_the_base_url_is_defined_once_in_the_header_every_page_includes(self):
        header = read_template("header_template.html")
        self.assertIn("const API_BASE_URL = window.location.origin + '/api/v1';", header)
        self.assertIn("const ANNOTATIONS_API_URL = `${API_BASE_URL}/annotations`;", header)
        all_templates = ""
        for name in sorted(os.listdir(CODE_DIRECTORY)):
            if not name.endswith(".html"):
                continue
            template = read_template(name)
            all_templates += template
            if name.endswith("_template.html") and name != "header_template.html":
                self.assertIn('{% include "header_template.html" %}', template, name)
                # The QC pages request their data through qc_shared_js.html's drawQCPageOnLoad
                self.assertTrue("API_BASE_URL" in template or "drawQCPageOnLoad(" in template, name)
        self.assertIn("const response = await fetch(`${API_BASE_URL}/${endpoint}`);",
                      js_function_body(read_template("qc_shared_js.html"), "function drawQCPageOnLoad("))
        self.assertEqual(1, all_templates.count("window.location.origin + '/api/v1'"))

    def test_the_results_page_reads_annotations_from_the_shared_url(self):
        # refreshRowAnnotations and fetchLocusTags, both through fetchJsonOrNull, which logs a
        # failed request and returns null
        self.assertEqual(2, read_template("index_page_template.html").count(
            "await fetchJsonOrNull(`${ANNOTATIONS_API_URL}/${encodeURIComponent(locusId)}`, "))


class GeneRegionListTests(unittest.TestCase):
    """Both pages take the gene regions from one list, header_template.html's GENE_REGIONS, and draw
    all four checkbox blocks (include and exclude, on each page) from one macro."""

    def test_the_list_matches_the_server_and_the_shared_checkboxes(self):
        regions = js_string_list(read_template("header_template.html"), "GENE_REGIONS")
        server = read_template("results_server.py")
        server_map = server[server.index("GENE_REGION_MAP = {"):]
        self.assertEqual(regions, re.findall(r'^    "(\w+)": \[', server_map[:server_map.index("\n}\n")], re.MULTILINE))
        checkboxes = filter_control_macros().gene_region_checkboxes
        self.assertEqual(regions, re.findall(r'name="exclude_gene_region_(\w+)"', str(checkboxes("exclude_gene_region_"))))

    def test_all_four_checkbox_blocks_render_the_macro(self):
        for name, prefix in (("index_page_template.html", ""), ("swim_plot_template.html", "swim_")):
            template = read_template(name)
            self.assertIn('{%- from "filter_control_macros.html" import gene_region_checkboxes, population_filter_field %}',
                          template, name)
            for kind in ("gene_region_", "exclude_gene_region_"):
                self.assertEqual(1, template.count("{{ gene_region_checkboxes('%s%s') | indent(" % (prefix, kind)),
                                 (name, kind))

    def test_both_pages_loop_over_it(self):
        # applyStateToUI and readStateFromUI, each once per include and exclude prefix of
        # GENE_REGION_FILTERS; getSwimParams once each for include and exclude
        template = read_template("index_page_template.html")
        self.assertIn("const GENE_REGION_FILTERS = [['gene_regions', 'gene_region_'], "
                      "['exclude_gene_regions', 'exclude_gene_region_']];", template)
        self.assertIn("for (const r of GENE_REGIONS) {", js_function_body(template, "function applyStateToUI("))
        self.assertIn("currentState[key] = GENE_REGIONS\n            .filter(",
                      js_function_body(template, "function readStateFromUI("))
        self.assertEqual(2, read_template("swim_plot_template.html").count("GENE_REGIONS.forEach(r => {"))


class FullTextPopupTests(unittest.TestCase):
    """One function writes the click-to-see-the-full-text span the .motif-truncated handler reads."""

    def test_every_shortened_motif_and_locus_id_goes_through_it(self):
        template = read_template("index_page_template.html")
        self.assertEqual(1, template.count('class="motif-truncated"'))
        self.assertIn('class="motif-truncated"', js_function_body(template, "function renderFullTextPopupTrigger("))
        row = js_function_body(template, "function renderLocusRow(locus, rowIndex) {")
        self.assertIn("formatMotifForDisplay(motifMain, motifTruncLen)", row)
        self.assertIn("formatMotifForDisplay(canonicalMotif, motifTruncLen)", row)
        self.assertIn("renderFullTextPopupTrigger(locusIdRaw, locusIdText, 'Click to see full locus ID')", row)


class DetailPanelMotifTests(unittest.TestCase):
    """Motifs over 30 bp are shortened in the detail panel's sample rows and affected-allele list."""

    def test_both_places_use_the_shortening_helper(self):
        detail = js_function_body(read_template("index_page_template.html"),
                                  "function renderLocusDetail(detail, outlierType) {")
        self.assertEqual(2, detail.count("formatMotifForDisplay(locus.Motif)"))
        self.assertNotIn("escapeHtml(locus.Motif || '')", detail)


if __name__ == "__main__":
    unittest.main()
