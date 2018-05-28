var chart;

/**
 * Request data from the server, add it to the graph and set a timeout
 * to request again
 */
function requestData() {
    $.ajax({
        url: '/live-data',
        success: function (points) {
            console.log(points)
            var series = chart.series[0];
            //var shift = series.data.length > 20; // shift if the series is longer than 20
            var shift = false;
            console.log(points.length)
            // add the point
            //for (var i = 0; i < points.length; i++) {
            //    point = points[i]
            //    //console.log(point)
            //    chart.series[0].addPoint(point, true, shift);
            //    
            //}
            console.log(chart.series[0].data)
            // call it again after one second
            setTimeout(requestData, 1000);
        },
        cache: false
    });
}

$(document).ready(function () {
    chart = new Highcharts.Chart({
        chart: {
            renderTo: 'data-container',
            defaultSeriesType: 'spline',
            events: {
                load: requestData
            }
        },
        title: {
            text: 'Temperature process'
        },
        xAxis: {
            type: 'datetime',
            //dateTimeLabelFormats: { // don't display the dummy year
            //    month: '%e. %b',
            //    year: '%b'
            //},
            title: {
                text: 'Datetime'
            }
        },
        yAxis: {
            title: {
                text: 'Temperature [C]'
            }
        },

        legend: {
            layout: 'vertical',
            align: 'right',
            verticalAlign: 'middle',
            enabled: false
        },

        plotOptions: {
            series: {
                label: {
                    connectorAllowed: false
                },
            }
        },

        series: [{
            name: 'Temperature',
            data: temp
        }],

        responsive: {
            rules: [{
                condition: {
                    maxWidth: 500
                },
                chartOptions: {
                    legend: {
                        layout: 'horizontal',
                        align: 'center',
                        verticalAlign: 'bottom'
                    }
                }
            }]
        }
    });
});