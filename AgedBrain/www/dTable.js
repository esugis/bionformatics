$(document).ready(function(){
    $(".shiny-datatable-output").on("mouseenter", "tbody", function() {
        tag();
    }).on("mouseleave", "tbody", function() {
        $(".tag").unbind() ;
    });
});

this.tag = function() {
    var my_t = "" ;
    $(".tag").mouseenter(function(event){
        //alert("entered here!") ;
            my_t = this.title ;
            this.title = "" ;
        var offset = $(this).offset() ;
        var tarr = my_t.split(",") ;
        var rarr = $.map(tarr, function(el,i) {
            var tag = el.split(":") ;
            return("<span style=\"font-size:" + ((Math.log(parseFloat(tag[1]))/-35) * 100 + 75) + "%\">" + tag[0] + "</span>")
        }) ;
        $("body").append("<div id=\"tooltip\"></div>") ;
        $("#tooltip").css({
            "top":offset.top + 40,
            "left":offset.left - 100
        }).html(rarr.join(" ")) ;
    }) ;
    $(".tag").mouseleave(function(){
        $("#tooltip").remove() ;
        this.title = my_t ;
    });
}

