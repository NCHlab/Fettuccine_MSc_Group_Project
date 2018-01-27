function imageclick() {
alert("Don't Click the image please.....");
}

  function format(value) {
      return '<div>Hidden Value: ' + value + '</div>';
  }
  $(document).ready(function () {
      var table = $('#Families').DataTable({});

      // Add event listener for opening and closing details
      $('#Families').on('click', 'td.details-control', function () {
          var tr = $(this).closest('tr');
          var row = table.row(tr);

          if (row.child.isShown()) {
              // This row is already open - close it
              row.child.hide();
              tr.removeClass('shown');
          } else {
              // Open this row
              row.child(format(tr.data('child-value'))).show();
              tr.addClass('shown');
          }
      });
  });
