
Highcharts.chart('press_chart', {

    title: {
        text: 'Pressure Chart'
    },

    subtitle: {
        text: 'Sensor 1'
    },

    yAxis: {
        title: {
            text: 'Units (S)'
        }
    },
    xAxis: {
        type: 'datetime',
        dateTimeLabelFormats: { // don't display the dummy year
            month: '%e. %b',
            year: '%b'
        },
        title: {
            text: 'Date'
        }
    },
    legend: {
        layout: 'vertical',
        align: 'right',
        verticalAlign: 'middle'
    },

    plotOptions: {
        series: {
            label: {
                connectorAllowed: false
            },

        }
    },

    series: [{
        name: 'Installation',
        data:
            { }
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