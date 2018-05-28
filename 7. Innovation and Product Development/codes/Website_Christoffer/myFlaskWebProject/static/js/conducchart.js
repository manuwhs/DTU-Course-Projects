Highcharts.chart('conduc_chart', {

    title: {
        text: 'Conductivity process'
    },

    //subtitle: {
    //    text: 'Sensor 1'
    //},

    yAxis: {
        title: {
            text: 'Conductivity'
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
        name: "Conductivity",
        data: conduc_data,
        color: '#FF9655'
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