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
               
	h3("Integrated data collection"),
	br(),
	p("AgedBrainSYSBIO database integrates public and private data resources as well as knowledge available for internal project specific analysis and collaborations.
        Integrated data set provides a snapshot of the current knowledge about the interactions in the Alzheimer's disease
         and gives the opportunity to study Alzheimer's disease through the network of complex interactions aggregated in one data source"),
	br(),
	p("Various data sources are integrated together in the form of a network. The connections in the networks correspond to the type of interaction.
	For each interaction thecode of origin describes the origin of the dataset. 
	The original interactors identificators(gene ID, proteins ID, SNP ID) are mapped to ENSEMBL gene ID."),
	p("Integrated data collection contains the following data types:"),
	br(),
	h3("Aggregated co-expression score"),
	br(),
	img(src='e_mem.png', align="left", width="70%",align="middle"),
	br(),br(),
	h3("Epistatic effects"), 
	br(),
	img(src='e_epi.png', align="left", width="70%",align="middle"),
	br(),br(),
	h3("Protein-protein interaction from IntAct database"),
	br(),
        img(src='e_intact.png', align="left", width="70%",align="middle"),
	br(),br(),
	h3("Protein-protein interactions from Hybrigenics"),
	br(),
	img(src='e_hyb.png', align="left", width="70%",align="middle"),
	br(),br(),
	h3("Summary of the interactions."),
	br(),
	img(src='edges.png', align="left", width="70%",align="middle"),
	p("Where :"), 
	p("ENSG1, ENSG2 - Ensemble gene identificator,"),
	p("PPI - protein-protein interaction,"),
	p("PCI - protein- protein complex interaction,"),
	p("IGRI - intergenic regions' interaction,"),
	p("HEI - Hybrigenics Experimental Interactions"),
	p("HLC - Hybrigenics Literature Collection"),
	p("IAHS - IntAct Homo Sapience"),
	p("ADIA - Alzheimer's Disease IntAct"),
	p("PDIA - Parkinson's Disease IntAct"),
	p("SIA - Synapse IntAct"),
	p("ADNI_VER - Alzheimer's Disease Neuroimaging Initiative and phenotype: ventricle enlargement rate,"), 
	p("TGEN - Translational Genomics Research Institute,"),
	p("HBTRC- Harvard Brain Tissue Resource Center,"),
	p("AND - Alzheimer's and non-Alzheimer's disease.
	"),
	br(),br(),
	h3("Node attributes"),
	br(),
	p("Each node is decribed by the set of attributes and the corresponding values: ENSG IDs, GWAS (p-value), 
	positive selection (p-value), gene type and the name of 221 brain tissues from Allen brain 
	atlas (Z-scores over mean expression values in one tissue are computed for each gene in the integrated dataset)."),
	img(src='nodes.png', align="left", width="70%",align="middle"),
	br(),br(),
	h3("BioModels"),
	br(),
	p("Collection of mathematical models from the literature (curated and non-curated) which are related to: Alzheimer's disease , 
	Huntington's disease , Parkinson's disease , amyotrophic lateral sclerosis , 
	neurodegenerative disease and tauopathy."),
	br(),br(),
	h3("Protein summaries from data resources in EMBL-EBI Proteomics Services Team"),
	br(),
	p("This data source provides aggregated information about selected proteins based on resources provided by EMBL-EBI Proteomics Services Team"),
	br(),br(),
	h3("Protein information"),
	br(),
	p("Protein sequence and functional information from Uniprot."),
	br(),br(),
	h3("Gene expression"),
	br(),
	p("1794 gene expression data sets deposited in Multi Experiment Matrix.",a("List of Datasets",
        href="http://biit.cs.ut.ee/mem/datacollections.cgi?db=mem_240212#A-AFFY-44",target="_blank")),

        p("24 Alzheimer specific publicly available gene expression data sets.",a("List of
        Datasets",href="http://biit.cs.ut.ee/mem/datacollections.cgi?db=agedbrain",target="_blank")),
	br(),br(),
		
	h3("GeneMania (release 24-Oct-2014, public resource)."),
	br(),
	p("GeneMANIA is a web tool that for identifying other genes that are related to a set of input genes, using a large set of functional association data.
	 Association data include protein and genetic interactions, pathways, co-expression, 
	co-localization and protein domain similarity."),
	br(),br(),		
	h4("REACTOME (v51, public resource).")
		
		),
		tabPanel('Help',

               		dataTableOutput("help"))#,
   )

  )#,
)
))


