/*jslint browser: true*/
/*global $*/


$(document).ready(function () {
  $('#serverside_table').DataTable({
    bProcessing: true,
    bServerSide: true,
    scrollX: true,
    sPaginationType: "full_numbers",
    bjQueryUI: true,
    sAjaxSource: '/tables/l1_serverside_table',
    columns: [
      {"data": "Column A"},
      {"data": "Column B"},
      {"data": "Column C"},
      {"data": "Column D"}
    ]
  });
});
