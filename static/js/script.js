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

    // Phylogenetic tree
    function genPhyloTree(data, info) {
        var basicNames = info.seq_names.basic_dataset;

        // Create tree
        var width = 800;
        var height = 2000;

        var svg = d3.select('#phylotree .scroll-wrapper')
            .append('svg')
            .attr('width', width)
            .attr('height', height)
            .append('g');

        // Create the cluster layout
        var cluster = d3.cluster()
            .size([height - 20, width - 200]);

        // Give the data to this cluster layout
        var root = d3.hierarchy(data, function (d) {
            return d.children;
        });
        cluster(root);

        // Add the links between nodes
        svg.selectAll('path')
            .data(root.descendants().slice(1))
            .enter()
            .append('path')
            .attr('d', function (d) {
                return 'M' + d.y + ',' + d.x
                    + 'L' + d.parent.y + ',' + d.x
                    + ' ' + d.parent.y + ',' + d.parent.x;
            })
            .style('fill', 'none')
            .attr('stroke', '#ccc');

        // Create nodes
        var g = svg.selectAll('g')
            .data(root.descendants())
            .enter()
            .append('g')
            .attr('transform', function (d) {
                return 'translate(' + d.y + ',' + d.x + ')'
            }).filter(function(d) { return d.data.name.indexOf('|') === -1; });

        // Add circles
        g.append('circle')
            .attr('r', 7);

        // Add text
        g.append('text')
            .attr('dy', '0.37em')
            .attr('x', 10)
            .text(function(d) { return d.data.name; })
            .attr('class', function(d) { return basicNames.indexOf(d.data.name) !== -1 ? 'basic' : 'related'; });
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
        }),
        $.getJSON('static/data/clustalo_phylotree.json', function(d) {
            data.tree = d;
        })
    ).then(function() {
        $('#msa .scroll-wrapper').append(genMSATable(data.alignment, data.info));
        genPhyloTree(data.tree, data.info);
    });

});
