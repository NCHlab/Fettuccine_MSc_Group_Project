/*jslint browser: true*/
/*global $*/

function jsconv(vars) {
	return vars
}

let myChart1 = document.getElementById('myChart_counts').getContext('2d');

// Global Options
Chart.defaults.global.defaultFontSize = 15;
Chart.defaults.global.defaultFontColor = '#777';

let barChart1 = new Chart(myChart1, {
	type:'doughnut', // bar, horizontalBar, pie, line, doughnut, radar, polarArea
	data:{
		labels:['ERV1', 'ERV1?', 'ERVK', 'ERVL', 'ERV-MaLR', 'ERVL?'],
		datasets:[{
			label:'Population',
			data:[
				jsconv({{arr_sc[0]}}),
				jsconv({{arr_sc[1]}}),
				jsconv({{arr_sc[2]}}),
				jsconv({{arr_sc[3]}}),
				jsconv({{arr_sc[4]}}),
				jsconv({{arr_sc[5]}})
			],
			//backgroundColor:'green',
			backgroundColor:["#FAA43A", "#F15854", "#F17CB0", "#60BD68", "#5DA5DA", "#F4F11F"],
			borderWidth:1,
			borderColor:'#777',
			hoverBorderWidth:3,
			hoverBorderColor:'#000'
		}]
	},
	options:{
		title:{
			display:true,
			text:'Genome Frequency',
			fontSize:20
		},
		legend:{
			display:false,
			position:'bottom',
			labels:{
				fontColor:'#000',
			usePointStyle: true,
			horizontalAlign: "center",
			verticalAlign: "center"
			}
		},
		layout:{
			padding:{
				left:30,
				right:0,
				bottom:0,
				top:0
			}
		},
		tooltips:{
			enabled:true
		}
	}
});

let myChart2 = document.getElementById('myChart_prots').getContext('2d');

let barChart2 = new Chart(myChart2, {
	type:'doughnut', // bar, horizontalBar, pie, line, doughnut, radar, polarArea
	data:{
		labels:['ERV1', 'ERV1?', 'ERVK', 'ERVL', 'ERV-MaLR', 'ERVL?'],
		datasets:[{
			label:'Population',
			data:[
				jsconv({{arr_sc_2[0]}}),
				jsconv({{arr_sc_2[1]}}),
				jsconv({{arr_sc_2[2]}}),
				jsconv({{arr_sc_2[3]}}),
				jsconv({{arr_sc_2[4]}}),
				jsconv({{arr_sc_2[5]}})
			],
			//backgroundColor:'green',
			backgroundColor:["#FAA43A", "#F15854", "#F17CB0", "#60BD68", "#5DA5DA", "#F4F11F"],
			borderWidth:1,
			borderColor:'#777',
			hoverBorderWidth:3,
			hoverBorderColor:'#000'
		}]
	},
	options:{
		title:{
			display:true,
			text:'Protein Frequency',
			fontSize:20
		},
		legend:{
			display:false,
			position:'bottom',
			labels:{
				fontColor:'#000',
			usePointStyle: true,
			horizontalAlign: "center",
			verticalAlign: "center"
			}
		},
		layout:{
			padding:{
				left:30,
				right:0,
				bottom:0,
				top:0
			}
		},
		tooltips:{
			enabled:true
		}
	}
});
