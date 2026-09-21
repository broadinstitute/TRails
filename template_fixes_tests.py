"""Regression tests for the HTML template findings (W15-W18, X13-X17, Y6, Y19-Y22, Z4, Z12-Z15,
V1, V13-V15) from the code review.

The templates are Jinja + JavaScript, so there is no Python harness that can execute them.
These tests pin the template source text instead, which is enough to catch a silent revert of
each fix.
"""

import os
import re
import unittest

REPO_DIRECTORY = os.path.dirname(os.path.abspath(__file__))


def read_template(name):
    """Returns the text of a file in the repository root.

    Args:
        name: File name, e.g. "index_page_template.html" or "results_server.py".

    Returns:
        The file's contents as a string.
    """
    with open(os.path.join(REPO_DIRECTORY, name)) as handle:
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


def js_handler_body(text, signature):
    """Returns the source of a top-level jQuery handler registration, closing `});` included.

    Args:
        text: Template source text to search.
        signature: Start of the registration, e.g. "$(document).on('click', '.add-note-btn".

    Returns:
        The text from the registration through the `});` in column zero that ends it.
    """
    start = text.index(signature)
    end = text.index("\n});\n", start)
    return text[start:end + 5]


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


class TagFilterRestoreTests(unittest.TestCase):
    """Z4: the saved tag must be selected in a dropdown whose menu already exists.

    applyStateToUI selects the tag by value, and Semantic UI silently ignores a 'set selected'
    for a value the menu does not hold yet, so the Search and Export handlers (which re-read the
    dropdown through readStateFromUI) would drop the tag on the next request.
    """

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_startup_waits_for_both_dropdowns_before_restoring_url_state(self):
        self.assertIn("Promise.all([tagDropdownReady, populateSampleIdDropdown()]).then(() => {\n"
                      "        // Read URL state and perform initial search after dropdowns are populated\n"
                      "        readStateFromUrl();", self.template)
        # The restore must not be gated on the sample-id list alone.
        self.assertNotIn("populateSampleIdDropdown().then(", self.template)

    def test_the_tag_population_promise_is_kept_rather_than_fired_and_forgotten(self):
        self.assertIn("const tagDropdownReady = populateTagDropdown();", self.template)
        self.assertLess(self.template.index("const tagDropdownReady = populateTagDropdown();"),
                        self.template.index("Promise.all([tagDropdownReady,"))

    def test_search_and_export_still_read_the_tag_off_the_dropdown(self):
        # This is what makes the ordering matter, so pin it: if these ever stopped reading the
        # dropdown the fix above would no longer be load bearing.
        self.assertIn("currentState.tag = $('#tag-filter').dropdown('get value') || '';", self.template)
        for handler in ("$('#search-button').on('click', function() {", "async function performExport(format) {"):
            self.assertIn("readStateFromUI();", self.template[self.template.index(handler):][:2000])


class SortIconTests(unittest.TestCase):
    """Z12: every sortable column header needs the span the active-sort styling colors."""

    def test_all_sortable_headers_carry_a_sort_icon(self):
        template = read_template("index_page_template.html")
        headers = re.findall(r"<th[^>]*\bdata-sort=\"([^\"]+)\"[^>]*>(.*?)</th>", template)
        self.assertEqual({"region", "size", "count", "variation_cluster"},
                         {field for field, _ in headers})
        for field, body in headers:
            self.assertIn('<span class="sort-icon">', body, field)

    def test_the_active_sort_class_is_only_styled_through_that_span(self):
        # Why the span is required: header_template.html colors .sort-icon, nothing else.
        self.assertIn("#results-table thead th.sorted-active .sort-icon", read_template("header_template.html"))


class ActiveSortClassTests(unittest.TestCase):
    """W15: the sort indicator class must not claim a direction the server does not apply.

    A header click only adds or removes the field from the sort list; each field's direction is
    fixed in results_server.build_api_order_by's SORT_MAPPING, and Region sorts ascending. The
    old markup styled `sorted-desc` for every sorted column and left `sorted-asc` unmatched.
    """

    def test_no_template_uses_a_direction_named_sort_class(self):
        for name in sorted(os.listdir(REPO_DIRECTORY)):
            if name.endswith(".html"):
                template = read_template(name)
                self.assertNotIn("sorted-asc", template, name)
                self.assertNotIn("sorted-desc", template, name)

    def test_the_class_the_indicator_applies_is_the_class_the_stylesheet_matches(self):
        indicator = js_function_body(read_template("index_page_template.html"),
                                     "function syncSortHeaderIndicator(")
        self.assertIn("removeClass('sorted-active')", indicator)
        self.assertIn("addClass('sorted-active')", indicator)
        styled = re.findall(r"#results-table thead th\.([\w-]+) \.sort-icon",
                            read_template("header_template.html"))
        self.assertEqual(["sorted-active"], styled)

    def test_the_server_still_sorts_the_region_column_ascending(self):
        # This is what made a single "sorted-desc" class wrong for every sorted column.
        self.assertIn('"region": ("gene_region_rank", "ASC")', read_template("results_server.py"))


class DetailRowCellSpecificityTests(unittest.TestCase):
    """W16: the detail-row cell rule must not reach the tables nested inside that cell.

    `#results-table tbody tr.locus-detail-row td` (1 id, 1 class, 3 types) out-specified
    `#results-table .sample-subtable td` (1,1,1), `.threshold-marker` (0,1,0) and
    `.invisible-table td` (0,1,1), so every sample cell got the outer cell's 40px indent,
    #fafafa background and 2px bottom border.
    """

    def setUp(self):
        self.header = read_template("header_template.html")

    def test_the_detail_cell_rule_uses_a_child_combinator(self):
        self.assertIn("#results-table tbody tr.locus-detail-row > td {", self.header)
        self.assertNotIn("#results-table tbody tr.locus-detail-row td {", self.header)

    def test_the_detail_row_still_has_exactly_one_direct_cell(self):
        # The child combinator is only equivalent for the outer cell because the row is built
        # with a single colspan td that the detail HTML is dropped into.
        template = read_template("index_page_template.html")
        self.assertIn('<tr class="locus-detail-row"><td colspan="18"', template)
        self.assertEqual(1, template.count('<tr class="locus-detail-row">'))

    def test_the_invisible_table_rule_outranks_the_main_table_cell_rule(self):
        # `.invisible-table td` alone loses to `#results-table tbody td`, which the Affected
        # Unsolved summary table sits inside.
        self.assertIn(".invisible-table td,\n#results-table .invisible-table td {", self.header)
        self.assertIn("#results-table tbody td {", self.header)

    def test_the_compact_sub_table_rules_are_still_the_ones_that_should_win(self):
        self.assertIn("#results-table .sample-subtable td { padding: 2px 6px;", self.header)


class AnnotationErrorCommentTests(unittest.TestCase):
    """W17 / V13: the comment must name the hook that actually returns the read-only-server 403.

    The comment is copied into both pages that write annotations, and W17 corrected only the
    results page, so V13 found the swim page still naming a function that does not exist.
    """

    TEMPLATES_WITH_THE_COMMENT = ("index_page_template.html", "swim_plot_template.html")

    def test_the_comment_names_an_existing_server_function(self):
        for name in self.TEMPLATES_WITH_THE_COMMENT:
            template = read_template(name)
            self.assertIn("(the results_server.restrict_writes before_request hook)", template, name)
            self.assertNotIn("require_local_writes", template, name)

    def test_no_template_names_a_server_symbol_that_does_not_exist(self):
        # Every `results_server.<name>` a comment cites has to be findable in the server.
        server = read_template("results_server.py")
        for name in self.TEMPLATES_WITH_THE_COMMENT:
            for symbol in set(re.findall(r"results_server\.([a-zA-Z_][a-zA-Z_0-9]*)", read_template(name))):
                self.assertTrue(re.search(r"^(def %s\(|%s = )" % (symbol, symbol), server, re.MULTILINE),
                                "%s cites results_server.%s, which does not exist" % (name, symbol))

    def test_that_function_exists_and_returns_the_detail_field_the_helper_reads(self):
        server = read_template("results_server.py")
        self.assertIn("def restrict_writes():", server)
        hook = server[server.index("def restrict_writes():"):]
        hook = hook[:hook.index("\n@app.")]
        self.assertIn('"detail": "This server is bound to a non-loopback address', hook)
        self.assertIn("403", hook)


class QC2OutlierWarningTests(unittest.TestCase):
    """W18: the per-source outlier map is gone from the QC2 page, which could never read it."""

    def test_the_page_no_longer_serializes_the_map(self):
        template = read_template("qc2_template.html")
        self.assertNotIn("OUTLIER_WARNINGS_BY_SOURCE", template)
        self.assertNotIn("outlier_warnings_by_source", template)

    def test_the_only_warning_call_on_the_page_passes_no_source(self):
        # `renderMendelianWarning()` with an empty argument list is the prose mention in the
        # comment above MENDELIAN_WARNINGS, not a call.
        calls = re.findall(r"renderMendelianWarning\(([^)]+)\)", read_template("qc2_template.html"))
        self.assertEqual(["item.sample_id"], calls)

    def test_the_shared_helper_tolerates_the_missing_const(self):
        # Deleting the const is only safe because the helper short-circuits on a falsy source and
        # guards the lookup with typeof.
        self.assertIn("const outlierMap = (source && typeof OUTLIER_WARNINGS_BY_SOURCE !== 'undefined')",
                      read_template("qc_shared_js.html"))

    def test_the_pages_that_do_pass_a_source_still_define_the_map(self):
        for name in ("index_page_template.html", "sample_qc_template.html", "swim_plot_template.html"):
            self.assertIn("const OUTLIER_WARNINGS_BY_SOURCE =", read_template(name), name)


class DeadSearchResultsGlobalTests(unittest.TestCase):
    """Z13: performSearch must not stash the response payload in a global nothing reads."""

    def test_no_template_writes_a_search_results_global(self):
        for name in sorted(os.listdir(REPO_DIRECTORY)):
            if name.endswith(".html"):
                self.assertNotIn("searchResults", read_template(name), name)


class AnnotationErrorBodyTests(unittest.TestCase):
    """Z14: annotation write failures must surface the server's own explanation."""

    def test_index_write_helpers_read_the_response_body(self):
        template = read_template("index_page_template.html")
        for signature in ("async function saveNote(", "async function deleteNote(",
                          "async function addTag(", "async function removeTag("):
            body = js_function_body(template, signature)
            self.assertIn("throw await annotationError(response);", body, signature)
            self.assertNotIn("new Error(`HTTP ${response.status}`)", body, signature)
        helper = js_function_body(template, "async function annotationError(")
        self.assertIn("await response.json().catch(() => ({}))", helper)
        self.assertIn("body.detail || body.error || `HTTP ${response.status}`", helper)

    def test_swim_plot_write_helpers_read_the_response_body(self):
        template = read_template("swim_plot_template.html")
        for signature in ("async function swimSaveNote(", "async function swimDeleteNote(",
                          "async function swimAddTag(", "async function swimRemoveTag("):
            body = js_function_body(template, signature)
            self.assertIn("throw await swimAnnotationError(response);", body, signature)
            self.assertNotIn("new Error('HTTP ' + response.status)", body, signature)
        helper = js_function_body(template, "async function swimAnnotationError(")
        self.assertIn("await response.json().catch(() => ({}))", helper)
        self.assertIn("body.detail || body.error || `HTTP ${response.status}`", helper)

    def test_the_server_still_sends_the_fields_the_helpers_read(self):
        server = read_template("results_server.py")
        self.assertIn('"error": "annotation changes are limited to this machine"', server)
        self.assertIn('"detail": "This server is bound to a non-loopback address', server)
        self.assertIn('{"error": "note_text is required and cannot be empty"}', server)


class ClipboardFallbackTests(unittest.TestCase):
    """Z15: the copy button must survive an unavailable or refused clipboard API."""

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
    """X13: the expand/collapse handler must not fire for the table's header cell."""

    def test_handler_is_bound_to_table_cells_only(self):
        template = read_template("index_page_template.html")
        self.assertIn("$(document).on('click', 'td.col-expand'", template)
        self.assertNotIn("$(document).on('click', '.col-expand'", template)

    def test_header_cell_still_carries_the_class(self):
        # The th keeps the class for column sizing, which is why the handler has to be scoped.
        self.assertIn('<th class="col-expand">', read_template("index_page_template.html"))


class TagAttributeReadTests(unittest.TestCase):
    """X16/X17: a tag named "true" or "1" must survive being read back off the badge.

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


class SampleQCSourceSelectorTests(unittest.TestCase):
    """X14: the unreachable database-source selector and its machinery are gone."""

    def test_selector_and_its_machinery_are_removed(self):
        template = read_template("sample_qc_template.html")
        for removed in ("qc-source-select", "qcRequestSeq", "sourceSelect", "db_labels"):
            self.assertNotIn(removed, template)

    def test_data_is_fetched_for_the_page_default_source(self):
        template = read_template("sample_qc_template.html")
        self.assertIn("sample_qc_data?source=${encodeURIComponent(DEFAULT_SOURCE)}", template)


class SharedJsHeaderTests(unittest.TestCase):
    """X15: qc_shared_js.html's header comment must match what the file defines and who uses it."""

    def setUp(self):
        self.shared_js = read_template("qc_shared_js.html")
        self.header = self.shared_js.split("=========================================================================")[1]

    def test_header_names_every_including_template(self):
        for name in sorted(os.listdir(REPO_DIRECTORY)):
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
    """Y6: the 'Unknown' motif-size bin the server emits must reach the charts and the popup."""

    def test_the_server_still_emits_an_unknown_bin(self):
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


class SampleQCChartCleanupTests(unittest.TestCase):
    """Y21: the chart-destruction bookkeeping is gone now that nothing re-renders a container."""

    def test_section_charts_bookkeeping_is_removed(self):
        template = read_template("sample_qc_template.html")
        self.assertNotIn("sectionCharts", template)
        self.assertNotIn("destroy()", template)


class SourceSelectorRemovalTests(unittest.TestCase):
    """Y22: the multi-database source selector is gone from the index and swim-plot pages."""

    def test_no_template_still_renders_a_source_selector(self):
        for name in ("index_page_template.html", "swim_plot_template.html", "sample_qc_template.html"):
            template = read_template(name)
            for removed in ("db_labels", 'name="source"', "swim_source", "available_sources"):
                self.assertNotIn(removed, template, name)

    def test_index_state_functions_no_longer_read_or_write_the_source_radio(self):
        template = read_template("index_page_template.html")
        self.assertNotIn("'require_above_population', 'source'", template)
        self.assertNotIn("currentState.source = ", template)
        # The default source is still what the API calls and the warning lookups use.
        self.assertIn("params.set('source', currentState.source);", template)

    def test_swim_plot_loads_the_default_source(self):
        self.assertIn("const params = { outlier_type: outlierType, source: SWIM_DEFAULT_SOURCE };",
                      read_template("swim_plot_template.html"))


class PopulationDataFilterGatingTests(unittest.TestCase):
    """Y19: include_loci_without_population_data only counts when require_above_population is set."""

    def test_both_endpoints_gate_the_flag(self):
        template = read_template("index_page_template.html")
        sends = list(re.finditer(r"params\.set\('include_loci_without_population_data'", template))
        # One send to /loci (buildApiUrl) and one to /export (performExport), each under a
        # require_above_population test in the statement that encloses it.
        self.assertEqual(2, len(sends))
        for send in sends:
            # The enclosing statement starts after the previous statement or block ends.
            statement = template[max(template.rfind(";", 0, send.start()),
                                     template.rfind("}", 0, send.start())):send.start()]
            self.assertIn("require_above_population", statement)

    def test_the_filter_chip_is_suppressed_with_no_population_threshold(self):
        self.assertIn("if ((key === 'population_metric' || key === 'include_loci_without_population_data')\n"
                      "            && !currentState.require_above_population) continue;",
                      read_template("index_page_template.html"))


class VariationClusterHelpTextTests(unittest.TestCase):
    """Y20: the size-diff help popup must not contradict the checkbox help below it."""

    def test_missing_cluster_is_described_as_unknown_not_as_absent(self):
        template = read_template("index_page_template.html")
        self.assertNotIn("have no variation beyond the repeat itself", template)
        self.assertIn("Loci with no variation cluster size are kept whatever value you enter, "
                      "because there is no size to compare against", template)


class SwimNoteEditingGuardTests(unittest.TestCase):
    """V1: the swim plot's Edit Note button must stay off until the stored note is known.

    showDotDetail opens the modal over a loading indicator and only then fetches the locus
    detail, and the footer handler takes the note's current text out of the rendered content.
    A click while the fetch was in flight therefore opened an empty editor over an existing
    note, and saving replaced that note with whatever was typed.
    """

    def setUp(self):
        self.template = read_template("swim_plot_template.html")
        self.show_dot_detail = js_function_body(self.template, "function showDotDetail(item) {")

    def test_the_footer_button_starts_disabled(self):
        # No locus is loaded before the first dot is clicked, so the markup itself is disabled.
        self.assertIn('<button class="ui basic button disabled swim-edit-note-btn" '
                      'style="float: left;" disabled>', self.template)

    def test_editing_is_disabled_before_the_modal_opens(self):
        self.assertLess(self.show_dot_detail.index("setSwimNoteEditingEnabled(false)"),
                        self.show_dot_detail.index("$('#dot-detail-modal').modal('show')"))

    def test_editing_is_re_enabled_only_after_the_content_is_rendered(self):
        self.assertEqual(1, self.show_dot_detail.count("setSwimNoteEditingEnabled(true)"))
        self.assertLess(self.show_dot_detail.index("renderDotDetailContent(item, detail, stats, params)"),
                        self.show_dot_detail.index("setSwimNoteEditingEnabled(true)"))

    def test_a_stale_response_cannot_re_enable_editing(self):
        # The newer-dot guard returns before the enable, so a late response for a dot the user
        # has already navigated away from cannot unlock the button for the one now loading.
        self.assertLess(self.show_dot_detail.index("if (seq !== swimDetailSeq) return;"),
                        self.show_dot_detail.index("setSwimNoteEditingEnabled(true)"))

    def test_a_failed_load_leaves_editing_disabled(self):
        catch_block = self.show_dot_detail[self.show_dot_detail.index(".catch(err =>"):]
        self.assertIn("Failed to load details: ", catch_block)
        self.assertNotIn("setSwimNoteEditingEnabled(true)", catch_block)

    def test_the_helper_disables_the_element_itself_not_only_its_styling(self):
        helper = js_function_body(self.template, "function setSwimNoteEditingEnabled(enabled) {")
        self.assertIn("swimDetailLoaded = enabled;", helper)
        self.assertIn(".toggleClass('disabled', !enabled).prop('disabled', !enabled)", helper)

    def test_the_click_handler_refuses_while_the_note_is_still_unknown(self):
        handler = js_handler_body(self.template, "$(document).on('click', '.swim-edit-note-btn'")
        self.assertIn("if (!_currentModalItem || !swimDetailLoaded) return;", handler)
        # The guard has to precede the read of the rendered note and the save.
        self.assertLess(handler.index("!swimDetailLoaded"), handler.index(".note-text"))


class ResultsPageNoteLoadFailureTests(unittest.TestCase):
    """V1: the results page must not open an empty note editor when the load failed.

    The results page fetches the note before opening the modal, so it has no click-while-loading
    window, but it used to treat both a non-ok response and a thrown request as "no note" and
    open an empty textarea. Saving that would have replaced the stored note the same way.
    """

    def setUp(self):
        self.handler = js_handler_body(read_template("index_page_template.html"),
                                       "$(document).on('click', '.add-note-btn, .note-indicator'")

    def test_a_failed_response_is_not_read_as_an_absent_note(self):
        self.assertNotIn("response.ok ? await response.json() : {}", self.handler)
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
    """V14: the detail panel's variation-cluster line must show the same three states as the row.

    A locus whose variation cluster was computed and then filtered out carries a filter reason
    and a NULL size difference, so gating the whole line on the size difference made the panel
    read as "no variation cluster" for exactly the loci whose VC cell shows a crossed-out icon.
    """

    def setUp(self):
        template = read_template("index_page_template.html")
        self.row = js_function_body(template, "function renderLocusRow(locus, rowIndex) {")
        self.detail = js_function_body(template, "function renderLocusDetail(detail, outlierType) {")

    def test_the_row_cell_still_has_both_states(self):
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

    def test_the_unfiltered_state_is_still_a_size_in_bp(self):
        self.assertIn("Variation Cluster: ${vcsSign}${vcsDiff.toLocaleString('en-US')}bp", self.detail)


class DetailRowOuterCellTests(unittest.TestCase):
    """V15: the detail row's own cell, never the hundreds of cells nested inside it.

    The rendered detail holds a per-sample sub-table plus the threshold, phenotype-score and
    affected-unsolved tables, so `.find('td')` on the row matched every one of those cells.
    refreshDetailPanel assigned the whole detail string to each, and the collapse handlers
    wrapped and animated each, firing their slideUp callback once per element.
    """

    def setUp(self):
        self.template = read_template("index_page_template.html")

    def test_no_detail_row_lookup_uses_a_descendant_selector(self):
        # Comment lines are allowed to name the old selector; code lines are not.
        for number, line in enumerate(self.template.split("\n"), start=1):
            if line.strip().startswith("//"):
                continue
            self.assertNotIn("find('td')", line, "line %d still selects descendant cells" % number)
            self.assertNotIn("$td.find('.slide-wrapper')", line,
                             "line %d still selects descendant wrappers" % number)

    def test_every_detail_cell_lookup_targets_a_direct_child(self):
        # Toggle off, the loading row's two renders, collapse one, collapse all, and the
        # post-annotation refresh.
        self.assertEqual(6, self.template.count(".children('td')"))
        # One slide wrapper per row: toggle off, collapse one, collapse all.
        self.assertEqual(3, self.template.count(".children('.slide-wrapper')"))

    def test_the_wrapped_cell_is_the_one_that_is_animated(self):
        self.assertIn("$td.wrapInner('<div class=\"slide-wrapper\"></div>');\n"
                      "    $td.children('.slide-wrapper').slideUp(200, function() {", self.template)

    def test_the_detail_row_still_has_exactly_one_direct_cell(self):
        # .children('td') is only equivalent to the intended target because the row is built
        # with a single colspan cell that the detail HTML is dropped into.
        self.assertIn('<tr class="locus-detail-row"><td colspan="18"', self.template)
        self.assertEqual(1, self.template.count('<tr class="locus-detail-row">'))


if __name__ == "__main__":
    unittest.main()
