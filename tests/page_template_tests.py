"""Tests that pin the source text of the results server's page templates.

The templates are Jinja + JavaScript, so there is no Python harness that can execute them.
These tests pin the template source text instead, which is enough to catch a silent change to
the behavior each test describes.
"""

import os
import re
import unittest

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
        for declaration in ("const tagDropdownReady = refreshTagFilterDropdown($('#tag-filter'));",
                            "const sampleIdDropdownReady = fetchJsonOrNull(`${API_BASE_URL}/schema`, 'Load schema')"):
            self.assertIn(declaration, self.template)
            self.assertLess(self.template.index(declaration), self.template.index("Promise.all([tagDropdownReady,"))

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
            for signature in self.HELPER_SIGNATURES + ("async function annotationError(",):
                self.assertIn(signature, self.header)
                self.assertNotIn(signature, template, name)
            for call in ("openNoteEditor(", "addTag(", "removeTag("):
                self.assertIn(call, template, name)
        self.assertIn("saveNote(currentLocusId, noteText)", js_function_body(self.header, "async function openNoteEditor("))

    def test_the_write_helpers_read_the_response_body(self):
        for signature in self.HELPER_SIGNATURES:
            body = js_function_body(self.header, signature)
            self.assertIn("throw await annotationError(response);", body, signature)
            self.assertNotIn("new Error(`HTTP ${response.status}`)", body, signature)
        helper = js_function_body(self.header, "async function annotationError(")
        self.assertIn("await response.json().catch(() => ({}))", helper)
        self.assertIn("body.detail || body.error || `HTTP ${response.status}`", helper)

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

    def test_the_helper_reads_the_outlier_map_only_without_mendelian_only(self):
        # Leaving the map out of the QC2 page is only safe because the helper does not read it then
        helper = js_function_body(read_template("qc_shared_js.html"), "function renderMendelianWarning(")
        self.assertIn("function renderMendelianWarning(sampleId, { mendelianOnly = false } = {}) {", helper)
        self.assertIn("const outlier = mendelianOnly ? undefined : OUTLIER_WARNINGS[sampleId];", helper)


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
        self.assertIn("bins.forEach(bin => { binData[bin] = []; });", template)
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
        load = js_function_body(template, "async function loadQC2Data(")
        self.assertEqual(2, load.count("_BINS, TRIO_SECTION);"))
        self.assertEqual(2, load.count("_BINS, CENTER_SECTION);"))

    def test_the_center_section_titles_name_the_server_trio_minimum(self):
        # The cutoff is the server's MIN_TRIOS_PER_CENTER, which the qc2() route passes in
        template = read_template("qc2_template.html")
        self.assertEqual(2, template.count("centers with ≥{{ min_trios_per_center }} trios</div>"))
        self.assertIn("min_trios_per_center=MIN_TRIOS_PER_CENTER,", read_template("results_server.py"))

    def test_the_shared_builder_hands_the_box_plot_plugin_its_statistics(self):
        builder = js_function_body(read_template("qc_shared_js.html"), "function createSingleBinBoxPlotChart(")
        self.assertIn("chart.config._boxStats = computeBoxStats(values);", builder)
        self.assertIn("p95: computePercentile(sortedValues, 95),", builder)
        self.assertIn("p99: computePercentile(sortedValues, 99)", builder)


class PopulationDataFilterGatingTests(unittest.TestCase):
    """include_loci_without_population_data only counts when require_above_population is set."""

    def test_the_results_query_gates_the_flag(self):
        template = read_template("index_page_template.html")
        sends = list(re.finditer(r"params\.set\('include_loci_without_population_data'", template))
        # One send, in buildApiUrl, whose query both /loci and the locus export use, under a
        # require_above_population test in the statement that encloses it.
        self.assertEqual(1, len(sends))
        self.assertIn("params.set('include_loci_without_population_data'",
                      js_function_body(template, "function buildApiUrl("))
        for send in sends:
            # The enclosing statement starts after the previous statement or block ends.
            statement = template[max(template.rfind(";", 0, send.start()),
                                     template.rfind("}", 0, send.start())):send.start()]
            self.assertIn("require_above_population", statement)

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


class SharedNoteEditorTests(unittest.TestCase):
    """Both pages edit notes in the one multi-line note editor in header_template.html.

    A textarea keeps a note's line breaks, which a window.prompt text field cannot hold.
    """

    def setUp(self):
        self.header = read_template("header_template.html")

    def test_the_editor_is_a_textarea_in_the_shared_header(self):
        self.assertIn('<div id="edit-note-modal" class="ui small modal">', self.header)
        self.assertIn('<textarea id="note-textarea" rows="4"', self.header)
        for name in ("index_page_template.html", "swim_plot_template.html"):
            self.assertNotIn('id="edit-note-modal"', read_template(name), name)

    def test_both_pages_open_it(self):
        self.assertIn("openNoteEditor(locusId, refreshAfterNoteChange);", read_template("index_page_template.html"))
        swim = read_template("swim_plot_template.html")
        self.assertIn("openNoteEditor($(this).data('locus-id'), refreshSwimModal);", swim)
        self.assertNotIn("prompt('Enter note:'", swim)

    def test_it_opens_over_the_swim_plot_detail_modal(self):
        self.assertIn("allowMultiple: true,", js_function_body(self.header, "async function openNoteEditor("))

    def test_every_open_sets_the_editor_and_no_close_resets_it(self):
        # A reset when the hide animation ends could land after a quick reopen had already filled
        # the editor with the stored note, leaving it empty over that note.
        editor = js_function_body(self.header, "async function openNoteEditor(")
        self.assertIn("$textarea.val(data.note.note_text);", editor)
        self.assertNotIn("onHidden: function", editor)

    def test_save_and_delete_tell_the_page_to_repaint(self):
        editor = js_function_body(self.header, "async function openNoteEditor(")
        self.assertIn("noteEditorOnNoteChanged = onNoteChanged;", editor)
        self.assertIn("noteEditorOnNoteChanged(currentLocusId);", editor)
        delete_handler = self.header[self.header.index("$(document).on('click', '#delete-note-btn'"):]
        self.assertIn("noteEditorOnNoteChanged(locusId);", delete_handler[:delete_handler.index("\n});\n")])


class NoteEditorLoadFailureTests(unittest.TestCase):
    """The note editor must not open empty when the stored note could not be loaded.

    It fetches the note before opening the modal, so it has no click-while-loading window, but
    treating a non-ok response or a thrown request as "no note" would open an empty textarea, and
    saving that would replace the stored note with whatever was typed.
    """

    def setUp(self):
        self.handler = js_function_body(read_template("header_template.html"), "async function openNoteEditor(")

    def test_a_failed_response_is_not_read_as_an_absent_note(self):
        self.assertIn("if (!response.ok) throw new Error(`HTTP ${response.status}`);", self.handler)

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

    def test_the_row_cell_has_both_states(self):
        self.assertIn("} else if (locus.VariationClusterFilterReason) {", self.row)
        self.assertIn("times circle outline icon", self.row)

    def test_the_detail_panel_has_the_same_second_state(self):
        self.assertIn("} else if (locus.VariationClusterFilterReason) {", self.detail)
        self.assertIn("times circle outline icon", self.detail)

    def test_the_detail_panel_reuses_the_row_cell_wording_and_colors(self):
        for name in ("VC_FILTER_REASON_TITLES", "VC_FILTER_REASON_COLORS"):
            self.assertIn(name, self.row)
            self.assertIn(name, self.detail)

    def test_both_guard_the_size_difference_the_same_way(self):
        guard = "locus.VariationClusterSizeDiff != null && locus.VariationClusterSizeDiff !== ''"
        self.assertIn(guard, self.row)
        self.assertIn(guard, self.detail)
        self.assertIn("Number.isFinite(", self.detail)

    def test_the_unfiltered_state_is_a_size_in_bp(self):
        self.assertIn("Variation Cluster: ${vcsSign}${vcsDiff.toLocaleString('en-US')}bp", self.detail)


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
        # Toggle off, the loading row's error message, the detail render (loadLocusDetailInto,
        # shared by the expand handler and the post-annotation refresh), collapse one and collapse
        # all, plus the two colspans (a locus detail row's and a grouped row's) that count the
        # cells of the row they open under.
        self.assertEqual(7, self.template.count(".children('td')"))
        self.assertEqual(2, self.template.count("colspan=\"${$row.children('td').length}\""))
        # One slide wrapper per row: toggle off, collapse one, collapse all.
        self.assertEqual(3, self.template.count(".children('.slide-wrapper')"))

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

    def test_no_colspan_is_a_hard_coded_number(self):
        self.assertEqual([], re.findall(r'colspan="\d+"', self.template))

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
        self.assertIn('"family_count": (f"NumAffectedUnsolvedFamiliesAboveUnaffected_{ot}", "DESC")',
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
        for signature in ("function applyStateToUI(", "function readStateFromUI(", "function buildApiUrl("):
            self.assertIn("'min_repeats_above_population'", js_function_body(self.template, signature), signature)
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

    def test_the_server_accepts_min_pli_for_api_callers(self):
        self.assertIn('request.args.get("min_pli")', read_template("results_server.py"))

    def test_a_url_with_a_parameter_the_page_does_not_use_is_normalized_in_place(self):
        # readStateFromUrl keeps only DEFAULT_STATE keys, and then always rewrites the arrival
        # hash in place, so a URL carrying a parameter the page does not use (e.g. min_pli)
        # neither keeps it nor traps Back.
        read_state = js_function_body(self.template, "function readStateFromUrl(")
        self.assertIn("if (key in DEFAULT_STATE) {", read_state)
        self.assertIn("    replaceStateInUrl();\n\n    applyStateToUI();", read_state)

    def test_the_metric_stays_in_the_url_because_the_column_uses_it(self):
        state_hash = js_function_body(self.template, "function buildStateHash(")
        self.assertNotIn("key === 'population_metric'", state_hash)


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
        self.assertIn('<input type="radio" name="require_above_population" value="all">', self.template)
        self.assertIn('<input type="radio" name="swim_require_above_population" value="all">',
                      read_template("swim_plot_template.html"))
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
        self.assertIn("if (!response.ok) throw new Error(`HTTP ${response.status}`);", loader)
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
        for name, template in (("index_page_template.html", index), ("swim_plot_template.html", swim)):
            self.assertNotIn("['TenK10K', 'TenK10K']", template, name)


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
            for match in re.finditer(r"(\$\('#[\w-]+'\))\.dropdown\(\{(.*?)\n    \}\);", read_template(name), re.DOTALL):
                if "allowAdditions: true" in match.group(2):
                    dropdowns_with_additions.append(match.group(1))
                    self.assertIn("...ESCAPING_ADDITIONS_DROPDOWN_SETTINGS,", match.group(2), match.group(1))
        self.assertEqual(["$('#sample-id-filter')", "$('#add-tag-dropdown')", "$('#swim-sample-id-filter')"],
                         dropdowns_with_additions)

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


class SearchErrorResponseTests(unittest.TestCase):
    """A failed search whose response is not JSON (e.g. a proxy's error page) reports its HTTP status."""

    def test_the_error_body_is_read_with_a_fallback(self):
        search = js_function_body(read_template("index_page_template.html"), "async function performSearch(")
        self.assertIn("const errorData = await response.json().catch(() => ({}));\n"
                      "            throw new Error(errorData.detail || `HTTP ${response.status}`);", search)


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
                self.assertIn("API_BASE_URL", template, name)
        self.assertEqual(1, all_templates.count("window.location.origin + '/api/v1'"))

    def test_the_results_page_reads_annotations_from_the_shared_url(self):
        # refreshRowAnnotations and fetchLocusTags, both through fetchJsonOrNull, which logs a
        # failed request and returns null
        self.assertEqual(2, read_template("index_page_template.html").count(
            "await fetchJsonOrNull(`${ANNOTATIONS_API_URL}/${encodeURIComponent(locusId)}`, "))


class GeneRegionListTests(unittest.TestCase):
    """Both pages take the gene regions from one list, header_template.html's GENE_REGIONS."""

    def test_the_list_matches_the_server_and_both_pages_checkboxes(self):
        regions = js_string_list(read_template("header_template.html"), "GENE_REGIONS")
        server = read_template("results_server.py")
        server_map = server[server.index("GENE_REGION_MAP = {"):]
        self.assertEqual(regions, re.findall(r'^    "(\w+)": \[', server_map[:server_map.index("\n}\n")], re.MULTILINE))
        for name, prefix in (("index_page_template.html", ""), ("swim_plot_template.html", "swim_")):
            template = read_template(name)
            for kind in ("gene_region_", "exclude_gene_region_"):
                self.assertEqual(regions, re.findall(r'name="%s%s(\w+)"' % (prefix, kind), template), (name, kind))

    def test_both_pages_loop_over_it(self):
        # applyStateToUI and readStateFromUI, each for include and exclude; getSwimParams likewise
        self.assertEqual(4, read_template("index_page_template.html").count("for (const r of GENE_REGIONS) {"))
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
