var main = function () {
  var
    start = 0, end = 200,

    update_view = function(data, stat) {
      var text, i;
      // genomes
      if ( 'genomes' in data ) {
        text = '';
        for ( i in data['genomes'] ) {
          text += "<option value='" + data['genomes'][i]['location'] + "'>" + data['genomes'][i]['title'] + "</option>";
        }
        $('#genomes').html( "<select>" + text + "</select>" );
      }
      else {
        $('#genomes').html( "No genomes" );
      }
      // sequence
      if ( 'alignments' in data ) {
        text = '';
        for ( i in data['alignments'] ) {
          text += "<option value='" + data['alignments'][i]['location'] + "'>" + data['alignments'][i]['title'] + "</option>";
        }
        $('#alignments').html( "<select>" + text + "</select>" );
      }
      else {
        $('#alignments').html( "No alignments" );
      }
      // log
    },

    request_status = function() {
      $.ajax({
        type: 'POST',
        url: '/status',
        success: update_view,
        dataType: 'json'
      });
    },

    update_visualization = function(data, stat) {
      $('#reference').html( data['ref'] )
      $('#reference_position').html( start + " to " + end );
      $('#visualize_position').val(start);
    }

    visualize = function() {
      $.ajax({
        type: 'POST',
        url: '/visualize',
        data: { 'start': start, 'end': end, 'ref': $('#genomes select').val(), 'sam': $('#alignments select').val() },
        success: update_visualization,
        dataType: 'json'
      });
    },

    visualize_position = function() {
      start = $('#visualize_position').val();
      end = parseInt(start) + 200;
      visualize();
    };

  return {
    init: function() {
      request_status();
      $('#visualize_button').click( visualize );
      $('#visualize_position_button').click( visualize_position );
    }
  }
};
