library(shiny)
#library(DT)
#require(devtools)
#library(rCharts)
shinyServer(function(input, output) {

    # table mem_self with indication of net. of origin 
#    help= read.table(file="/group/inbox/nikolaeva/datatableshelp.txt", header=TRUE, sep="\t", quote="\"")
#    mem_self=read.table(file="/home/nikolaeva/AgedBrain/ppi2_data/TablePPISelfInt.displayIA.txt", header=TRUE, sep="\t")
#    mem_net=read.table(file="/home/nikolaeva/AgedBrain/ppi2_data/TablePPI.displayIA.txt",header=TRUE, sep="\t")
#    epi=read.table(file="/home/nikolaeva/AgedBrain/ppi2_data/epistatic_bin-bin_interactions_at_MAF_prod_GE_0.01_bonf_pvalue_LE_1.0E-11.tsv",header=TRUE,sep="\t")
#    proteins=read.table(file="/home/nikolaeva/AgedBrain/ppi2_data/proteinsUniprot.txt", header=TRUE,sep="\t")
load(file = "/home/nikolaeva/AgedBrain/ppi2_data/integrated_ds/int_with_nodes11022016.RData") 
load(file = "/home/nikolaeva/AgedBrain/ppi2_data/ABSBShinyData.RData") 
   # Lemps's data (self interacting PPIs) customize the length drop-down menu; display 10 rows per page by default

    output$table5 = renderDataTable({
        mem_self[, input$show_vars_self, drop = FALSE]
    },options = list(aLengthMenu = c(10, 30, 50, 100, 250), iDisplayLength = 10), escape = FALSE)

   # Added network customize the length drop-down menu; display 10 rows per page by default
    output$table6 = renderDataTable({
      mem_net[, input$show_vars_net, drop = FALSE]
    },options = list(aLengthMenu = c(10, 30, 50, 100, 250), iDisplayLength = 10), escape = FALSE)

   # Added epistatic interaction pairs customize the length drop-down menu; display 10 rows per page by default
   # output$table7 = renderDataTable({
   #     epi[, input$show_vars_epi, drop = FALSE]
   # },options = list(aLengthMenu = c(10, 30, 50, 100, 250), iDisplayLength = 10), escape = FALSE)

       # Added proteins info customize the length drop-down menu; display 10 rows per page by default
    output$table8 = renderDataTable({
        proteins
    },options = list(aLengthMenu = c(10, 30, 50, 100, 250), iDisplayLength = 10), escape = FALSE)

   # Added network customize the length drop-down menu; display 10 rows per page by default
    
	output$table9 = renderDataTable({
      	net[, input$show_vars_intnet, drop = FALSE]
    	},options = list(aLengthMenu = c(10, 30, 50, 100, 250), iDisplayLength = 10), escape = FALSE)

  



   # Added network customize the length drop-down menu; display 10 rows per page by default
    output$help = renderDataTable({
        help
    },escape = FALSE)
   #web pageopening
   getPage<-function(){
    return((HTML(readLines('http://www.google.com'))))
  }
  output$inc<-renderUI({
   # x <- input$ebi  
    getPage()
  })
datasetInput <- reactive({
    switch(input$dataset,
#	   "Integrated data collection" = net,
	   "MEM" = mem_net,
           "MEM self interactions" = mem_self,
           "Epistatic interactions" = epi,
	   "Integrated data collection" = net
           )
  })
output$downloadData <- downloadHandler(
    filename = function() {
                 paste(input$dataset, '.csv', sep='')
         },
    content = function(file) {
      write.csv(datasetInput(), file)
    }
  )


})
