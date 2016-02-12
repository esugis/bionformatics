$(document).ready(function(){
  $(".shiny-datatable-output").on("mouseenter","table tr", function(){
   alert("table");

});
});

/*$(function() {
  setTimeout(function() {
    Shiny.unbindAll();
    // Do your work and then call Shiny.bindAll()
  }, 10);
});*/



/*$(document).ready(function(){  
 Shiny.unbindAll();
$("#table1").children("div").mouseenter(function(){
Shiny.unbindAll();  
alert("Yes!Table!");
Shiny.bindAll();
 });

});*/  

/*$( window ).load(function() {
    alert($("#table1").length)

});*/
var refreshIntervalId = setInterval(function () {
/*Shiny.unbindAll();*/
/*$("#DataTables_Table_0").mouseenter=alert("yes");*/ 
if($("#DataTables_Table_0 tbody tr").length){
	clearInterval(refreshIntervalId);
Shiny.unbindAll();
$("#DataTables_Table_0 tbody").onmouseenter=alert("Yes!Table!");
/*Shiny.bindAll();*/  

}
/*$("#DataTables_Table_0").mouseenter(function(){
    alert("Yes!Table!");
  });*/

    
 console.log($("#DataTables_Table_0"));
   }, 50);
