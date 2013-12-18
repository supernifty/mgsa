var main = function () {
  var
    update_view = function(data, stat) {
      var text, i;
      // genomes
      if ( 'genomes' in data ) {
        text = '';
        for ( i in data['genomes'] ) {
          text += "<li>" + data['genomes'][i]['title'] + "</li>";
        }
        $('#genomes').html( "<ul>" + text + "</ul>" );
      }
      else {
        $('#genomes').html( "No genomes" );
      }
      // sequence
      if ( 'sequences' in data ) {
        text = '';
        for ( i in data['sequences'] ) {
          text += "<li>" + data['sequences'][i]['title'] + "</li>";
        }
        $('#sequences').html( "<ul>" + text + "</ul>" );
      }
      else {
        $('#sequences').html( "No sequence" );
      }
      // log
      if ( 'log' in data ) {
        text = '';
        for ( i in data['log'] ) {
          text += "<li>" + data['log'][i]['when'] + " " + data['log'][i]['what'] + "</li>";
        }
        $('#log').html( "<ul>" + text + "</ul>" );
      }
      else {
        $('#log').html( "No log entries" );
      }
    },

    request_status = function() {
      $.ajax({
        type: 'POST',
        url: '/status',
        success: update_view,
        dataType: 'json'
      });
    }
  return {
    init: function() {
      request_status();
    }
  }
};
