Highcharts.chart('temp_chart', {

    title: {
        text: 'Temperature process'
    },

    //subtitle: {
    //    text: 'Sensor 1'
    //},

    yAxis: {
        title: {
            text: 'Temperature [C]'
        }
    },
    xAxis: {
        type: 'datetime',
        dateTimeLabelFormats: { // don't display the dummy year
            month: '%e. %b',
            year: '%b'
        },
        title: {
            text: 'Datetime'
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
        name: "Temperature",
        data: temp_data
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