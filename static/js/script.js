$(function() {

    // Tab selection
    $('.button').on('click', function() {
        var sel = $(this).data('tab');
        $('.tab, .button').removeClass('active');
        $('#' + sel).addClass('active');
        $(this).addClass('active');
    });

    // Multiple sequence alignment table
    function genMSATable(tableData, tableInfo) {
        var basicNames = tableInfo.seq_names.basic_dataset,
            relatedNames = tableInfo.seq_names.related_dataset;

        var segmentStart = tableInfo.optimal_segment.start,
            segmentEnd = tableInfo.optimal_segment.end;

        var table = $('<table></table>');
        $(tableData).each(function(i, rowData) {
            var row = $('<tr></tr>');
            $(rowData).each(function(j, cellData) {
                // Modify row names
                if (cellData === 'markings') {
                    cellData = 'Clustal Omega Markings'
                } else if (cellData === 'num_unique') {
                    cellData = 'Number of Unique Amino Acids'
                }

                // Add colors
                var td = '<td';

                if (basicNames.indexOf(cellData) !== -1) td += ' class="basic"';
                if (relatedNames.indexOf(cellData) !== -1) td += ' class="related"';

                // AA color scheme based on Clustal Omega
                if (['A', 'V', 'F', 'P', 'M', 'I', 'L', 'W'].indexOf(cellData) !== -1) td += ' class="red"';
                if (['D', 'E'].indexOf(cellData) !== -1) td += ' class="blue"';
                if (['R', 'H', 'K'].indexOf(cellData) !== -1) td += ' class="magenta"';
                if (['S', 'T', 'Y', 'H', 'C', 'N', 'G', 'Q'].indexOf(cellData) !== -1) td += ' class="green"';

                // Highlight the segment of alignment that is retained (shift 1 col because row names)
                if (segmentStart < j && j <= segmentEnd) td += ' bgcolor="#eee8aa"';

                td += '>';

                row.append($(td + cellData + '</td>'));
            });
            table.append(row);
        });
        return table;
    }

    // Get data to display
    var data = {};
    $.when(
        $.ajax({
            type: 'GET',
            url: 'static/data/clustalo_alignment.csv',
            success: function(d) {
                data.alignment = Papa.parse(d.trim()).data;
            }
        }),
        $.getJSON('static/data/dataset_info.json', function(d) {
            data.info = d;
        })
    ).then(function() {
        $('#msa .scroll-wrapper').append(genMSATable(data.alignment, data.info));
    });

});
