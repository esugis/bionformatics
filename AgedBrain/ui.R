library(shiny)

#    help= read.table(file="/group/inbox/nikolaeva/datatableshelp.txt", header=TRUE, sep="\t", quote="\"")
#    mem_self=read.table(file="/home/nikolaeva/AgedBrain/ppi2_data/TablePPISelfInt.displayIA.txt", header=TRUE, sep="\t")
#    mem_net=read.table(file="/home/nikolaeva/AgedBrain/ppi2_data/TablePPI.displayIA.txt",header=TRUE, sep="\t")
#    epi=read.table(file="/home/nikolaeva/AgedBrain/ppi2_data/epistatic_bin-bin_interactions_at_MAF_prod_GE_0.01_bonf_pvalue_LE_1.0E-11.tsv",header=TRUE,sep="\t")
#    proteins=read.table(file="/home/nikolaeva/AgedBrain/ppi2_data/proteinsUniprot.txt",	header=TRUE,sep="\t")
#    load integrated dataset net.		
load(file = "/home/nikolaeva/AgedBrain/ppi2_data/integrated_ds/int_with_nodes11022016.RData")
load(file = "/home/nikolaeva/AgedBrain/ppi2_data/ABSBShinyData.RData")

shinyUI(fluidPage(
  headerPanel(
  'AgedBrainSYSBIO data collection',
        tags$head(
	tags$img(src='Logo_AgedBrainSysbio.png', align="left", width="18%"),
        tags$img(src='banniere.jpg',align="left"),
        tags$script(src='dTable.js'),
        tags$link(rel='stylesheet', type='text/css', href='dTable.css')
       )
  ),
 title = 'AgedBrainSYSBIO data collection',
  sidebarLayout(
    sidebarPanel(width=3,
	conditionalPanel(
        'input.dataset === "Introduction"',
	 p('Introduction')
     	 ),

     	conditionalPanel(
        'input.dataset === "Integrated data collection"',
#        checkboxGroupInput('show_vars_intnet', 'Choose columns to display:',
#                           names(net), selected = names(net))
 	 checkboxGroupInput('show_vars_intnet',p('Choose columns to display:'),
                           names(net), selected = c("ENSG.A","ENSG.B","SCORE","INTERACTION",
							"EVIDENCE.ORIGIN","GENE.NAME.A","P.VALUE.A",
							"GENE.NAME.B","P.VALUE.B"))
#,
#        downloadButton('downloadData_intnet', 'Download')
     	),
	
	conditionalPanel(
        'input.dataset === "MEM"',
        checkboxGroupInput('show_vars_net', p('Choose columns to display:'),
                           names(mem_net), selected = c("P1","P2","mem_p","tag","DS_p","mem_n","DS_n","GeneMANIA","IntAct"))
#                           names(mem_net), selected = names(mem_net))
#,
#	downloadButton('downloadData_net', 'Download')
      ),

        conditionalPanel(
        'input.dataset === "MEM self interacting PPI"',
        checkboxGroupInput('show_vars_self',p( 'Choose columns to display:'),
                           names(mem_self), selected = c("P1","P2","mem_p","tag","DS_p","mem_n","DS_n","GeneMANIA","IntAct"))
#                           names(mem_self), selected = names(mem_self))
#,
#	downloadButton('downloadData_self', 'Download')
      ),

#      conditionalPanel(
#        'input.dataset === "Epistatic interactions"',
#	checkboxGroupInput('show_vars_epi', 'Choose columns to display:',
#                          names(epi), selected = names(epi)),
#	downloadButton('downloadData_epi', 'Download')
#      ),

	conditionalPanel(
        'input.dataset === "UniProt"',
        p('Links to the protein information page in UniProt'),
	img(src='uniprot_logo_trans.png', align="left", width="50%",align="middle")
      ),
	conditionalPanel(
        'input.dataset === "Protein summaries&Pathways"',
        p('Links to the EMBL-EBI BioModels,protein summaries from data resources in EMBL-EBI  Proteomics Services Team
	 and the description of the selected pathways by The Babraham Institute'),
	img(src='EMBL_EBI_Logo_black.png', align="left", width="50%",align="middle"),
	img(src='BI-logo_3_1.png', align="left", width="50%",align="middle")
      ),
	conditionalPanel(
        'input.dataset === "Data"',
        p('List of data used in the analysis')
      ),
	conditionalPanel(
        'input.dataset === "Help"',
        p('Discribes meaning of the values in the columns of MEM and MEM self-interactions')
      )
    ),



 mainPanel(
    # title = 'AgedBrainSYSBIO data collection',
    tabsetPanel(
      id = 'dataset',
     
	tabPanel('Introduction',style = "overflow:hidden;", 
	h4("Welcome to AgedBrainSYSBIO database!"),
	br(),
	h4("Data collection"),p("AgedBrainSYSBIO database integrates public and private data resources as well as knowledge available for internal project specific analysis and collaborations.
	Integrated data set provides a snapshot of the current knowledge about the interactions in the Alzheimer's disease
	 and gives the opportunity to study Alzheimer's disease through the network of complex interactions aggregated in one data source."),
	br(),
	img(src='combine_ds.png', align="left", width="70%",align="middle"),
       br(),br(),
	h4("Data integration"),p("The usage of the extended AgedBrainSYSBIO data collection together with analysis tools provides the researchers with possibility 
	to identify most relevant protein-protein interactions and their interacting partners associated with neurodegenerative disorders 
	such as Alzheimer's disease."),
        br(),br(),
	img(src='mapt_module.png', align="left", width="70%",align="middle"),
	br(),br(),
	h4("Statistical analysis of the complex interaction network"),p("Statistical analysis of the complex interaction network that aggregates various experimental interactions specific to 
	Alzheimer's disease provides the insights of the specific disease triggering mechanism."),
	br(),br(),
	img(src='node_degree.png', align="left", width="70%",align="middle"),
	br(),br(),br(),br(),br(),br(),br(),
	p(img(src='european_flag_blueyellow_standard.jpg',height="10%", width="8%",align="right"),span("AgedBrainSYSBIO has received funding from the European Union's Seventh Framework 
	Programme ",style="float: right; vertical-align: bottom;padding-right:10px"),span(" technological development and 
        demonstration under grant agreement No 
	305299 ",style="float: right; vertical-align: bottom;padding-right:10px"))
#	img(src='european_flag_blueyellow_standard.jpg', width="10%",align="right")
	),
	  tabPanel('Integrated data collection',
               dataTableOutput("table9")),
        tabPanel('MEM',
               dataTableOutput("table6")),
        tabPanel('MEM self interacting PPI',
               dataTableOutput("table5")),
#       tabPanel('Epistatic interactions',
#               dataTableOutput("table7")),
	tabPanel('UniProt',
	       dataTableOutput("table8")),
	tabPanel('Protein summaries&Pathways',
        br(),
	h4("AgedBrainSYSBIO: models related to neurodegenerative diseases"),
	br(),
	img(src='biomodels.png', align="left", width="70%",align="middle"),
	br(),br(),
	p("Collection of mathematical models from the literature (curated and non-curated)
 	which are related to: Alzheimer's disease , Huntington's disease , Parkinson's disease , amyotrophic lateral sclerosis , 
	neurodegenerative disease and tauopathy."),
	p("Please visit the AgedBrainSYSBIO collection of models related to neurodegenerative disease by following the link below:"),
	a("BioModels Database",href="http://www.ebi.ac.uk/biomodels-main/agedbrain",target="_blank"),
	br(),br(),
	h4("Protein summaries from data resources in EMBL-EBI Proteomics Services Team"),
	br(),
	img(src='protein_summary.png', align="left", width="70%",align="middle"),
#	img(src='EMBL_EBI_Logo_black.png', align="left", width="20%"),	
	br(),br(),
#	p(""),
	a("BIN1",href="http://wwwdev.ebi.ac.uk/AgedBrainSYSBIO/index.php?protein=bin1",target="_blank"),
	br(),
	a("TAU",href="http://wwwdev.ebi.ac.uk/AgedBrainSYSBIO/index.php?protein=tau",target="_blank"),
	br(),
	a("FYN",href="http://wwwdev.ebi.ac.uk/AgedBrainSYSBIO/index.php?protein=fyn",target="_blank"),
	br(),
	a("APP",href="http://wwwdev.ebi.ac.uk/AgedBrainSYSBIO/index.php?protein=app",target="_blank"),
	br(),
	a("LRRK2",href="http://wwwdev.ebi.ac.uk/AgedBrainSYSBIO/index.php?protein=lrrk2",target="_blank"),
	br(),br(),
	h4("Selected AgedBrainSYSBIO Pathways"),
	p("Relevant signalling pathways and metabolic networks were selected based on the existing
	knowledge on Alzheimer aetiology and symptomatology."),
	br(),
	img(src='BI-logo_3_1.png', align="left", width="30%",align="middle"),
	br(),br(),
	a("Glutamate Signalling Pathway",
	href="http://lenoverelab.org/projects/AgedBrainSYSBIO/pathways/Glutamate", target="_blank"),
	br(),
	a("Amyloid fibril formation pathway",href="http://lenoverelab.org/projects/AgedBrainSYSBIO/pathways/Amyloid",target="_blank"),
	br(),
        a("PIP3/AKT pathway",href="http://lenoverelab.org/projects/AgedBrainSYSBIO/pathways/PIP3AKT",target="_blank"),
	br(),
        a("MAP kinase pathway",href="http://lenoverelab.org/projects/AgedBrainSYSBIO/pathways/MAPK",target="_blank"),
	br(),
        a("Endocytic pathway",href="http://lenoverelab.org/projects/AgedBrainSYSBIO/pathways/endocytosis",target="_blank"),
	br(),
        a("Immune System",href="http://lenoverelab.org/projects/AgedBrainSYSBIO/pathways/immunesystem",target="_blank"),
	br(),
        a("Apoptotic cleavage of cellular proteins pathway",href="http://lenoverelab.org/projects/AgedBrainSYSBIO/pathways/Apoptosis",target="_blank"),
	br(),
        a("Oxidative stress induced senescence pathway",href="http://lenoverelab.org/projects/AgedBrainSYSBIO/pathways/OxidativeStress",target="_blank"),
	br(),
        a("Cell junction organization pathway",href="http://lenoverelab.org/projects/AgedBrainSYSBIO/pathways/cellJunction",target="_blank"),
	br(),
        a("Gap junction trafficking and regulation pathways",href="http://lenoverelab.org/projects/AgedBrainSYSBIO/pathways/GapJunction",target="_blank"),
	br()
	),

      tabPanel('Data',
	p("- Hybrigenics provided 2 datasets of experimental Y2H protein-protein intercations and 
	interactions collected from public databases. The collection is availbale from PimRider for AgedBrain  partners.
	Experimental Y2H data was detected in four different cDNA libraries: Human Fetal Brain (908 PPI), Human Brain (558 PPI),
	Drosophila head (13 PPI) and Mouse Adult Brain (18 PPI).  Data sets of protein-protein interactions (PPIs) 
	containing information about PPIs specific to DISC1 (publication PMID:17043677) and APP proteins (produced by HYBRIGENICS)."),
               
	p("- 1794 gene expression data sets deposited in Multi Experiment Matrix (provided by Quretec).",a("List of Datasets",
	href="http://biit.cs.ut.ee/mem/datacollections.cgi?db=mem_240212#A-AFFY-44",target="_blank")), 
	
	p("- Positive and negative co-expression scores for each protein pair computed using MEM over 1794 
	datsets including 77 datasets associated with brain (provided by Quretec). "),
	
	p("- Common key words (tag cloud) for protein pairs with especially strong correlation (or 
	anticorrelation) in the subset of datasets. Such subsets of datasets share common key words that describe 
	certain biological state or condition (generated by Quretec)."),
	
	p("- 24 Alzheimer specific publicly available gene expression data sets (collected by Quretec).",a("List of
	Datasets",href="http://biit.cs.ut.ee/mem/datacollections.cgi?db=agedbrain",target="_blank")),
	
	p("- Protein-protein interactions associated with Alzheimer's disease (PMID:21163940)."),
	
	p("- Protein profiling by protein array(PMID:21826230),downloaded to Quretec server."),

       	p("- List of GWAS genes (provided by INSERM)."),

	p("- GWAS list of genes from publications (collected by Quretec)."),
	
	p("- List of 248 positively selected genes in human European population (based on 1000G data) 
	(provided by INSERM). "),
		
	p("- Epistatic effects of pairs of SNPs on ventricle volume, detected in ADNI patients,presented as 
	a list of epistatic interactions between genomic regions, including both genes and intergenic 
	regions. The list is comprised a total of 2'943 interactions between 1'709 distinct regions. Considering only 
	protein-coding genes,735 genes are linked through 567 interactions(provided by SIB)."),
		
	p("- IntAct database (public data source, provided by EBI)."),
		
	p("- BioGrid ( ver.3.2.104, public resource)."),
		
	p("- GeneMania (release 24-Oct-2014, public resource)."),
		
	p("- STRING (ver. 9.1, public resource)."),
	
	p("- REACTOME (v51, public resource)."),
		
	p("- KEGG ( release 69.0, public resource).")
		),
		tabPanel('Help',

               		dataTableOutput("help"))#,
   )

  )#,
)
))


