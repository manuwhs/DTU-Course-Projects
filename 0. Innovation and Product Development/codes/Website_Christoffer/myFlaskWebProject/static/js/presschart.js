Highcharts.chart('press_chart', {

    title: {
        text: 'Pressure process'
    },

    //subtitle: {
    //    text: 'Sensor 1'
    //},

    yAxis: {
        title: {
            text: 'Pressure [Pa]'
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
        name: "Pressure",
        data: press_data,
        color: '#DDDF00'
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