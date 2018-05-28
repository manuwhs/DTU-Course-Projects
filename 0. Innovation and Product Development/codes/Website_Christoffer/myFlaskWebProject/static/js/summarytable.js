//sample data

var printIcon = function (value, data, cell, row, options) { //plain text value
    //console.log(value.cell.value)
    link = "<img class='infoImage' src='/static/icons/" + value.cell.value + ".gif' height='20'>";
    console.log(link)
    //return "<img class='infoImage' src='/static/icons/pass.png' height='25' width='20'>";
    return link;
};

var linktourl = function (value, data, cell, row, options) {
    console.log(value)
    url="result/" + value.cell.value
    return url;
};

$("#summary-table").tabulator({
    height: 600, // set height of table to enable virtual DOM
    data: tabledata, //load initial data into table
    layout: "fitColumns", //fit columns to width of table (optional)
    columns: [ //Define Table Columns
        //{ title: "A", field: "name", formatter: link, "www.google.de"},
        { title: "ID", field: "name", sorter: "string", width: 200, formatter: "link", formatterParams:{url:linktourl} },
        { title: "Datetime", field: "date", sorter: "string" },
        //{ title: "Status", field: "status", sorter: "string" },
        { title: "User", field: "user", sorter: "string" },
        { title: "temp-setValue [degree]", field: "temp", sorter: "number" },
        { title: "ph-setValue [no]", field: "ph", sorter: "string" },
        { title: "con-setValue[unit]", field: "con", sorter: "string" },
        { title: "Status", field: "status", formatter: printIcon, variableHeight: true}
        
    ],
    rowClick: function (e, id, data, row) { //trigger an alert message when the row is clicked
        row.select(); //toggle row selected state
        //alert("Row " + id + " Clicked!!!!");
    }
});