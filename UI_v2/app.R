options(stringsAsFactors = F, scipen = 999) 
library(ape)
library(stringr)
library(plotly)
library(ggtree)
library(shinythemes)
library(caper)
library(shinyWidgets)
library(shinyjs)
library(phylolm)
library(shinycssloaders)
library(RColorBrewer)
set_color = brewer.pal(8, 'Paired')
set_color = append(set_color,c("#fdfd99","#e2cc1a"))


std <- function(x) sd(x)/sqrt(length(x))

lm_eqn <- function(m=lm(Y ~ X,data)){
  # paste("R2 =", format(summary(m)$r.squared, digits = 2) , "; p-value = ",formatC(summary(m)$coefficients[2,4], format = "e", digits = 2))
  paste("R2 =", round(summary(m)$r.squared, 2) , "; p-value = ",formatC(summary(m)$coefficients[2,4], format = "e", digits = 0))
}


GLS <- function(dataframe=shorebird){
  aic = 1000000
  dt = data.frame()
  for (model in c("LM","lambda","OUfixedRoot","OUrandomRoot","BM")){
    for (measurement_error in c(T,F)){
      if (model == "LM"){
        fit = lm(pgls_y~pgls_x, data = dataframe$data)
        measurement_error = NA
      } else if (model != "lambda"){
        fit <- phylolm(pgls_y~pgls_x, phy = dataframe$phy, data = dataframe$data, model = model,measurement_error=measurement_error)
      } else{ fit <- phylolm(pgls_y~pgls_x, phy = dataframe$phy, data = dataframe$data, model = model)
      measurement_error = NA}
      a = summary(fit)
      if (length(a$optpar)==0){a$optpar=NA}
      if (length(a$aic)==0){a$aic=NA
      a$logLik=NA
      a$optpar=NA
      a$sigma2=NA}
      
      dt = rbind(dt,data.frame(
        model,
        measurement_error,
        p_val_slope = a$coefficients[2,4],
        r.squared = a$r.squared,
        adj.r.squared = a$adj.r.squared,
        aic = a$aic,
        logLik = a$logLik,
        optpar = a$optpar,
        sigma2 = a$sigma2
      ))
      if ( !is.na(a$aic < aic) & a$aic < aic ){ best_fit_model = fit
      best_model = model
      aic = a$aic}
    }
  }
  dt = dt[!duplicated(dt$aic),]
  dt = dt[order(dt$aic),]
  return(list(dt,best_fit_model,best_model))
}


Clade_color = c("Other Invertebrates"="#f5b48a","Lepido Diptera"="red","Other Vertebrates"="#A6CEE3","Other Insecta"="#FF7F00",
                Nematoda="#B2DF8A",Teleostei="#1F78B4",Hymenoptera="#ba8e18",Aves="#5b5b5b",Mammalia="#66281A",Embryophyta="#33A02C"
)

table_phylo = read.delim("/home/fbenitiere/ICASDB_app/UI_v2/www/phylogenetic_tree/tree_description.tab")
phylogenetic_trees = paste("www/phylogenetic_tree/",table_phylo$name,sep="")
names(phylogenetic_trees) = table_phylo$description

data_by_species = data.frame(species="")
for (file in list.files("www/species_informations_tables",full.names = T)){
  dt = read.delim(file,header = T)
  data_by_species = merge(dt,data_by_species, by.x = "species", by.y = "species", all.x = TRUE, all.y = TRUE)
}


dt_species = read.delim("www/name_species.tab",row.names = "full_names",header=T)
listNomSpecies = tapply(dt_species$row.names,dt_species$clades,function(x)  str_replace_all(x,"_"," "))

data_by_species$clade.qual = factor(dt_species[data_by_species$species,]$clades, levels = c("Embryophyta","Lepido Diptera","Hymenoptera",
                                                                                            "Other Insecta","Nematoda","Other Invertebrates",
                                                                                            "Mammalia","Aves","Teleostei","Other Vertebrates"))


axisInter = read.delim("www/variables_information_tables/inter_axis.tab",sep="\t")
axisInter_quantitative = axisInter[axisInter$quantitative,]
axisInter_qualitative = axisInter[!axisInter$quantitative,]


axisInter_list_qualitative = tapply(axisInter_qualitative$display_label,axisInter_qualitative$group,list)
axisInter_list_quantitative = tapply(axisInter_quantitative$display_label,axisInter_quantitative$group,list)

axisIntra=list("Gene expression (FPKM)"="FPKMgene", # list des axes X et Y possible
               "Intron expression (FPKM of the gene)"="FPKMintron",
               "SVR per gene"="SVRgene",
               "NSVR per gene"="NSVRgene",
               "SVR per intron"="SVRintron",
               "NSVR per intron"="NSVRintron",
               "No introns per gene"="IntronPerGene",
               "Intron length on average per gene (bp)"="averageLength",
               "Intron length (bp)"="IntronLength",
               "N1 per intron"="N1",
               "N1+N2 per intron"="N1+N2"
)


# UI interface
ui <- shinyUI(fluidPage(style = "background-color:#a7c2da;height: 100%;font-family: Economica; font-size: 21px",tags$style(HTML( #couleur des onglets
  "
    .tabbable > .nav > li > a {background-color: #529bde; font-size: 23px; font-family: 'Economica' ; color:white}
    .tabbable > .nav > li[class=active]  > a {background-color: #136dc0; font-size: 25px; font-family: 'Economica' ; color:white}
    .my-class {font-size: 21px;}
    ")),
  setSliderColor(c("#529bde","#529bde","#529bde","#529bde","#529bde","#529bde"), c(1, 2, 3, 4, 5,6)),shinyjs::useShinyjs(),
  tabsetPanel(id = "tabs",
              tabPanel("Inter-species graphics",
                       column(12,
                              column(1,offset = 1,h2("Y axis :")),
                              column(10,selectInput("y_inter", "",choices = axisInter_list_quantitative)),
                              column(10,class = "well",offset = 0,withSpinner(type = 4,color ="#136dc0",plotlyOutput("plot_inter", width = "100%", height = "100%"))),
                              column(2,class="well", 
                                     fluidRow(column(4,
                                                     materialSwitch(inputId = "boxplot_inter", label = h3("Boxplot"),status="primary")),
                                              column(5,
                                                     prettyCheckboxGroup( "scale_inter", shape = "round",animation="pulse",bigger=T,
                                                                          status = "primary",
                                                                          fill = TRUE,h3("Scales"),choices = list("X log10","Y log10"),
                                                                          selected = list(""),inline = T))),
                                     prettyRadioButtons("busco_inter",shape = "round",animation="pulse",
                                                        status = "primary",bigger=T,
                                                        fill = TRUE,h3("Busco dataset"),
                                                        choices = list("Eukaryota"="eukaryota","Metazoa"="metazoa","Embryophyta"="embryophyta",
                                                                       "None"="None"),selected = "eukaryota",inline = T,),
                                     selectInput("shape_inter", h3("Shape"), choices = c("none",axisInter_list_qualitative) , selected = "none"),
                                     pickerInput(inputId = "clades_inter",label = h3("Select/deselect clades"), 
                                                 choices = levels(data_by_species$clade.qual),selected = levels(data_by_species$clade.qual), 
                                                 choicesOpt = list(
                                                   style = rep_len("font-size: 21px", length(levels(data_by_species$clade.qual)))
                                                 ),
                                                 options = list(style = "my-class",`actions-box` = TRUE,`selected-text-format` = "count > 3"
                                                 ), 
                                                 multiple = TRUE ),
                                     materialSwitch(inputId = "pgls_inter", label = h3("PGLS"),status="primary",inline = T),
                                     selectInput("tree_inter", h3("Tree for PGLS"), choices = phylogenetic_trees, selected = "1"),
                                     dropdown(
                                       tags$h2("Parameters"),
                                       pickerInput("svr_class",h3("Introns SVR classes"), choices = c( "All major introns" = "all",
                                                                                                       "SVR > 5%" = "high",
                                                                                                       "SVR <= 5%" = "low"
                                       )),
                                       
                                       sliderInput("coverage_inter",h3("Minimal coverage (reads/bp)"),min = 0, max = 1000, value = 200),
                                       fileInput("upload_data", h3("Choose data File"),
                                                 accept = c(
                                                   "text/csv/tab",
                                                   "text/comma-separated-values,text/plain",
                                                   ".tab")),
                                       
                                       style = "unite", icon = icon("gears"),
                                       status = "primary", 
                                     )
                              ),
                              column(1,offset = 3,h2("X axis :")),
                              column(5,selectInput("x_inter", "",choices = axisInter_list_quantitative))
                       )
                       
              ),
              
              ## INTRA SPECIES GRAPHIC
              tabPanel("Intra-species graphics",
                       column(12,
                              column(2,offset = 1,selectInput("species_selected_intra", h2("Species studied"),choices =  listNomSpecies , selected = "Drosophila melanogaster")),
                              column(2, style='padding:10px;',imageOutput('species_image_intra',height=5 , width=5), div(style = "height:120px;")),
                              column(3,selectInput("y_intra", h2("Y axis"),choices = axisIntra , selected = 1)),
                              column(3,selectInput("x_intra", h2("X axis or Histogram"),choices = axisIntra , selected = 1)),
                              column(9,class = "well",offset = 1,withSpinner(type = 4,color ="#136dc0", plotlyOutput("plot_intra",height = "100%",width="100%"))),
                              column(width =2,class="well",  
                                     
                                     fluidRow(column(4,
                                                     materialSwitch(inputId = "histogram_intra", label = h3("Histogram"),status="primary")),
                                              column(5,
                                                     prettyCheckboxGroup( "scale_intra", shape = "round",animation="pulse",bigger=T,
                                                                          status = "primary",
                                                                          fill = TRUE,h3("Scales"),choices = list("X log10","Y log10"),
                                                                          selected = list(""),inline = T))),
                                     prettyRadioButtons("busco_intra",shape = "round",animation="pulse",
                                                        status = "primary",bigger=T,
                                                        fill = TRUE,h3("Busco dataset"),
                                                        choices = list("Eukaryota"="eukaryota","Metazoa"="metazoa","Embryophyta"="embryophyta",
                                                                       "None"="None"),selected = "None",inline = T,),
                                     dropdown( 
                                       sliderInput("svr_range_intra",h3("SVR range of the introns studied"),min = 0, max = 0.5, value =  c(0,0.5)),
                                       sliderInput("bin_intra",h3("Proportion of N by points (%)"),min = 0, max = 100, value =  10),style = "unite", icon = icon("gears"),
                                       status = "primary", 
                                     ),
                                     downloadButton("download_fpkm", "Genes Expression", style = "color: #fff; background-color: #27ae60; border-color: #fff;font-size: 21px;"),
                                     downloadButton("download_busco_id", "Busco identification", style = "color: #fff; background-color: #27ae60; border-color: #fff;font-size: 21px;"),
                                     downloadButton("download_svr", "Introns Alternative Splicing", style = "color: #fff; background-color: #27ae60; border-color: #fff;font-size: 21px;")
                              )
                       )
              ),
              
              ### GENE STRUCTURE
              tabPanel("Gene structure",
                       column(10,offset=1,class = "well",fluidRow(
                         column(2,offset = 1,
                                selectInput("species_gene_struct", h4("Species studied"), 
                                            choices = listNomSpecies , selected = "Drosophila melanogaster")),
                         column(2, imageOutput('species_image',height=10,width=10),
                                div(style = "height:100px;")),
                         
                         column(2,
                                prettyRadioButtons("gene_list",shape = "round",animation="pulse",
                                                   status = "primary",bigger=T,
                                                   fill = TRUE,h4("ID choice"),
                                                   choices = c("gene id"="gene_id","metazoa busco id" = "busco_id_metazoa",
                                                               "eukaryota busco id" = "busco_id_eukaryota",
                                                               "embryophyta busco id" = "busco_id_embryophyta"),selected = "gene_id",inline = T,)),
                         
                         column(2, selectInput("studied_gene", h4("Gene studied"),choices="")),
                         column(2,  dropdown(   
                           sliderInput("sliderscale", h4("Bp Bins on x axis"),
                                       min = 0, max = 100000, value = 1000),
                           style = "unite", icon = icon("gears"),
                           status = "primary", 
                         ))))
                       ,column(10,offset = 1,class = "well",withSpinner(type = 4,color ="#136dc0",
                                                                        plotlyOutput("structureGene",width = "100%", height = "100%")))),
              
              ### PHYLOGENETIC TREE
              tabPanel("Phylogenetic tree",
                       fluidRow(
                         column(2,offset = 1,selectInput("select_tree", h3("Tree selected"),choices = phylogenetic_trees,
                                                         selected = phylogenetic_trees[1])),
                         column(2,selectInput("layout_tree", h3("Layout"),
                                              choices = c("roundrect","ellipse","circular","equal_angle","daylight")
                         )),
                         column(2,sliderInput("spacing", h3("Spacing"),min = 0, max = 1, value = 0.05
                         )),
                         column(2,sliderInput("tip_size", h3("Tips size"),min = 0, max = 10, value = 2
                         ))),
                       fluidRow(
                         column(10,offset = 1, class = "well",
                                column(6,
                                       withSpinner(type = 4,color ="#136dc0",plotOutput("plot_principal", height = 900,
                                                                                        brush = brushOpts(
                                                                                          id = "plot2_brush",
                                                                                          resetOnNew = TRUE
                                                                                        )
                                       )))
                                ,
                                column(6,h4("Zoom"),
                                       # dropdown(
                                       
                                       withSpinner(type = 4,color ="#136dc0",plotOutput("plot_zoomed", height = 600)),
                                       column(3,offset = 3,sliderInput("tip_size_zoom", h4("Tips size zoomed"),min = 0, max = 10, value = 2
                                       )),
                                       column(3,sliderInput("spacing_zoom", h4("Spacing zoomed"),min = 0, max = 1, value = 0.05
                                       ))
                                       # ,                                  tags$h3("Parameters"),
                                       
                                       # style = "unite", icon = icon("zoom-in", lib = "glyphicon"),
                                       # status = "primary")
                                )
                         )
                       )
              )
  )
)
)



##################### Server logic ##################### ##################### ##################### ##################### ##################### 
server <- function(input, output,session) {
  
  # Tout le bloc permet de montrer des images lors du passage de souris sur interSpecies
  addResourcePath(prefix = "imgResources", directoryPath = "www/species_images/")
  
  hover_event <- reactive({
    event_data(event = "plotly_hover", source = "hoverplotsource")
  })
  
  unhover_event <- reactive({
    event_data(event = "plotly_unhover", source = "hoverplotsource")
  })
  
  hoverplotlyProxy <- plotlyProxy("plot_inter", session)
  
  observeEvent(unhover_event(), {
    hoverplotlyProxy %>%
      plotlyProxyInvoke("relayout", list(images = list(NULL)))
  })
  
  observeEvent(input$upload_data,ignoreNULL = T, {
    inFile <- input$upload_data
    print(inFile)
    print(inFile$datapath)
    dt = read.delim(inFile$datapath,header = T)
    print(dt)
    data_by_species <<- merge(dt,data_by_species, by.x = "species", by.y = "species", all.x = TRUE, all.y = TRUE)
    print(colnames(data_by_species))
    
    axisInter <<- rbind(axisInter,data.frame(group = "Uploaded","display_label"=colnames(dt),name_label=colnames(dt),description=NA,quantitative=T))
    axisInter_quantitative <<- axisInter[axisInter$quantitative,]
    axisInter_list_quantitative <<- tapply(axisInter_quantitative$display_label,axisInter_quantitative$group,list)
    updateSelectizeInput(session, "y_inter", choices = axisInter_list_quantitative,options = list(maxOptions = 3000))
    updateSelectizeInput(session, "x_inter", choices = axisInter_list_quantitative,options = list(maxOptions = 3000))
  })
  
  observeEvent(hover_event(), {
    dt = unlist(str_split(hover_event()$customdata,"_;_"))
    
    if ("Y log10" %in% input$scale_inter  & "X log10" %in% input$scale_inter){
      hoverplotlyProxy %>%
        plotlyProxyInvoke("relayout", list(images = list(
          list(
            source = dt[5],
            xref = "x",
            yref = "y",
            x = hover_event()$x,
            y = hover_event()$y,
            sizex = (log10(as.numeric(dt[2])) - log10(as.numeric(dt[1])))/6,
            sizey = (log10(as.numeric(dt[4])) - log10(as.numeric(dt[3])))/6,
            opacity = 1
          )
        )))
    } else {
      hoverplotlyProxy %>%
        plotlyProxyInvoke("relayout", list(images = list(
          list(
            source = dt[5],
            xref = "x",
            yref = "y",
            x = hover_event()$x,
            y = hover_event()$y,
            sizex = (as.numeric(dt[2]) - as.numeric(dt[1]))/6,
            sizey = (as.numeric(dt[4]) - as.numeric(dt[3]))/6,
            opacity = 1
          )
        )))
    }
  })
  
  ###### Fin du bloc
  observeEvent(input$boxplot_inter,ignoreNULL = FALSE,{
    if (input$boxplot_inter){
      updateSelectizeInput(session, "x_inter", choices = axisInter_list_qualitative,options=list(maxOptions=3000))
      updatePrettyCheckboxGroup(session, 
                                prettyOptions = list(shape = "round",animation="pulse",bigger=T,
                                                     status = "primary",
                                                     fill = TRUE), "scale_inter", choices = list("Y log10"))
      updateMaterialSwitch(session, "pgls_inter", value = F)
      shinyjs::disable("pgls_inter")
    } else {    
      updateSelectizeInput(session, "x_inter", choices = axisInter_list_quantitative,options=list(maxOptions=3000))
      updatePrettyCheckboxGroup(session, 
                                prettyOptions = list(shape = "round",animation="pulse",bigger=T,
                                                     status = "primary",
                                                     fill = TRUE),"scale_inter", choices = list("X log10","Y log10"))
      shinyjs::enable("pgls_inter")
    }
  })
  
  
  observeEvent(input$pgls_inter, {
    if (!input$pgls_inter){
      shinyjs::disable("tree_inter")
    } else {
      shinyjs::enable("tree_inter")}
  })
  
  
  output$plot_inter <- renderPlotly({
    ylabel = axisInter[axisInter$display_label == input$y_inter,]$name_label 
    xlabel = axisInter[axisInter$display_label == input$x_inter,]$name_label
    
    if ( grepl("buscodataset_",xlabel)){ xlabel = paste(xlabel,input$busco_inter,".quant",sep="")}
    if ( grepl("buscodataset_",ylabel)){ ylabel = paste(ylabel,input$busco_inter,".quant",sep="")}
    
    if ( grepl("svr_class_",xlabel)){ xlabel = str_replace(xlabel,"all",input$svr_class)}
    if ( grepl("svr_class_",ylabel)){ ylabel = str_replace(ylabel,"all",input$svr_class)}
    print(ylabel)
    print(xlabel)
    if ( input$busco_inter != "None"){ 
      minimum_coverage = paste("median_coverage_exon.buscodataset_",input$busco_inter,".quant",sep="")
      data_by_species = data_by_species[data_by_species[minimum_coverage] > input$coverage_inter , ]
    }
    
    data_by_species = data_by_species[data_by_species$clade.qual %in% input$clades_inter,]
    
    data_by_species = data_by_species[!is.na(data_by_species[,xlabel]) & !is.na(data_by_species[,ylabel]),]
    if ( input$pgls_inter ){
      arbrePhylo = read.tree(input$tree_inter)
      data_by_species = data_by_species[data_by_species$species %in% arbrePhylo$tip.label,]
    }
    
    data_by_species$img_local <- paste("imgResources/",data_by_species$species,".png",sep="")
    if (input$shape_inter != "none"){
      data_by_species$shape = unlist(data_by_species[axisInter[axisInter$display_label == input$shape_inter,]$name_label ])
    } else {
      data_by_species$shape = ""
    }
    
    
    if ( !input$boxplot_inter  ){
      data_by_species$custom_data = paste(min(data_by_species[,xlabel]),max(data_by_species[,xlabel]),
                                          min(data_by_species[,ylabel]),max(data_by_species[,ylabel]),data_by_species$img_local,sep="_;_")
    } else {
      data_by_species$custom_data = "custom_data"
    }
    
    
    p = ggplot(
      data_by_species,aes_string(x=xlabel,y=ylabel, text="species",
                                 customdata="custom_data") )+ 
      scale_fill_manual(values=Clade_color) +
      xlab(input$x_inter) + 
      ylab(input$y_inter) + theme_bw() + theme(
        axis.title.x = element_text(color="black", size=25,family="economica"),
        axis.title.y = element_text(color="black", size=25, family="economica"),
        axis.text.y =  element_text(color="black", size=20, family="economica"),
        axis.text.x =  element_text(color="black", size=20, family="economica"),
        title =  element_text(color="black", size=15, family="economica"),
        legend.text =  element_text(color="black", size=20, family="economica")
      ) + geom_point(aes(fill=clade.qual,shape=shape),size=4,alpha=0.7)
    p
    
    if ( "Y log10" %in% input$scale_inter  ){
      lm_y = log10(data_by_species[,ylabel])
      p = p + scale_y_log10()
    } else {
      lm_y = data_by_species[,ylabel]
    }
    
    if ( "X log10" %in% input$scale_inter  ){ 
      lm_x = log10(data_by_species[,xlabel])
      p = p + scale_x_log10()
    } else {
      lm_x = data_by_species[,xlabel]
    }
    
    if ( input$boxplot_inter  ){
      p = p + geom_boxplot(alpha=.1,fill="grey") + 
        geom_point(aes(fill=clade.qual,shape=shape),size=3,alpha=0.7)+ theme_bw() + theme(
          axis.title.x = element_text(color="black",angle = 50, size=25,family="economica"),
          axis.title.y = element_text(color="black", size=25, family="economica"),
          axis.text.y =  element_text(color="black", size=20, family="economica"),
          axis.text.x =  element_text(color="black", size=25,angle = 50, family="economica"),
          title =  element_text(color="black", size=15, family="economica"),
          legend.text =  element_text(color="black", size=20, family="economica")
        ) + ggtitle(paste("N=",nrow(data_by_species))) + theme(legend.position='none')
    } else if ( !input$boxplot_inter ){
      if (input$pgls_inter){
        shorebird <- comparative.data(arbrePhylo, 
                                      data.frame(species=data_by_species$species,
                                                 pgls_x=lm_x,
                                                 pgls_y=lm_y), species, vcv=TRUE)
        
        gls = GLS(shorebird)
        print( gls[[1]] )
        print(gls[[2]])
        p = p + geom_abline(slope=1,intercept=0,alpha = .6) +
          ggtitle(paste("N=",nrow(data_by_species)," / LM:",lm_eqn(lm(lm_y ~ lm_x)),
                        " / PGLS:",lm_eqn(pgls(pgls_y~pgls_x,shorebird)),
                        "/ Best",gls[[3]],":",lm_eqn(gls[[2]])))
      } else { 
        p = p +geom_abline(slope=1,intercept=0,alpha = .6)+ ggtitle(paste("N=",nrow(data_by_species)," / LM:",lm_eqn(lm(lm_y~lm_x))))
      }
    }
    
    ggplotly( p, tooltip = c("text") , 
              height = 1000*.8, 
              hoverinfo = 'none',
              source = "hoverplotsource"
    ) %>%
      event_register('plotly_hover') %>%
      event_register('plotly_unhover')
    
  })
  
  
  observe({
    if (input$tabs == "Intra-species graphics"){
      species = input$species_selected_intra
      
      species = str_replace_all(species," ","_")
      species_genes = read.delim(paste("www/per_species_data/",species,"/by_gene_analysis.tab",sep="") , header=T , sep="\t",comment.char = "#")
      rownames(species_genes) = species_genes$gene_id
      species_intron = read.delim(paste("www/per_species_data/",species,"/by_intron_major_overlap.tab",sep=""))
      species_intron$median_fpkm = species_genes[species_intron$gene_id,]$median_fpkm
      species_intron$svr = species_intron$splice_variant_rate
      species_intron$nsvr = species_intron$nonsplice_variant_rate
      species_intron$sum_n1 = species_intron$n1
      species_intron$sum_n2 = species_intron$n2_spl3 + species_intron$n2_spl5
      species_intron$sum_n3 = species_intron$n3_spl3 + species_intron$n3_spl5
      species_intron$geneprop_tabid = species_intron$gene_id
      species_intron$length = abs(species_intron$splice5 - species_intron$splice3)
      species_intron <<- species_intron
    }
  })
  
  
  output$species_image_intra <- renderImage({
    secies_name <- str_replace(input$species_selected_intra," ","_")
    list(src=paste("www/species_images/",
                   secies_name,".png",sep=""))
  },deleteFile=FALSE)
  
  
  output$plot_intra <- renderPlotly({
    species = str_replace(input$species_selected_intra," ","_")
    color_species = Clade_color[dt_species[species,]$clades]
    
    species_intron = species_intron[ species_intron$into_cds == "True",]
    
    
    if ("None" != input$busco_intra){
      busco_gene = read.delim(paste("www/per_species_data/",species,"/busco_to_gene_id_",input$busco_intra,sep=""))
      
      species_intron = species_intron[species_intron$gene_id %in% busco_gene$gene_id, ]
    }
    
    
    species_intron = species_intron[!is.na(species_intron$svr) & 
                                      species_intron$svr >= input$svr_range_intra[1] & 
                                      species_intron$svr < input$svr_range_intra[2],]
    
    NSVRgene = tapply(species_intron$sum_n3,species_intron$geneprop_tabid,sum) #calcul NSVR par gene
    SVRgene = tapply(species_intron$sum_n2,species_intron$geneprop_tabid,sum) #calcul SVR par gene
    intronGene = tapply(species_intron$sum_n1,species_intron$geneprop_tabid,sum) # calcul N1 par gene
    average = tapply(species_intron$length,species_intron$geneprop_tabid,mean)
    FPKMgene = tapply(species_intron$median_fpkm,species_intron$geneprop_tabid,mean)
    
    SVRgene = (SVRgene/(SVRgene+intronGene)) #calcul du taux de SVR par gène
    
    NSVRgene = (NSVRgene/(NSVRgene + 2*intronGene)) #calcul du taux de NSVR par gène
    
    listeAxis = list("FPKMgene"=FPKMgene,"SVRgene"=SVRgene,"NSVRgene"=NSVRgene,"IntronPerGene"=table(species_intron$geneprop_tabid),"averageLength"=average,
                     "N1"=species_intron$sum_n1,"N1+N2"=species_intron$sum_n1+species_intron$sum_n2,
                     "SVRintron"=species_intron$svr,"FPKMintron"=species_intron$median_fpkm,"NSVRintron"=species_intron$nsvr,"IntronLength"=species_intron$length)
    
    if ( !input$histogram_intra ){
      xaxis = unlist(listeAxis[input$x_intra])
      proportion = input$bin_intra / 100
      quantile = unique(quantile(xaxis, probs = seq(0, 1,proportion),na.rm=T))
      intervalle = cut(xaxis, quantile,include.lowest = T,include.higher=T)
      
      X = tapply(xaxis, intervalle, mean)
      XerrorBar = tapply(xaxis, intervalle, std)
      yaxis = unlist(listeAxis[input$y_intra])
      Y = tapply(yaxis, intervalle, mean)
      YerrorBar = tapply(yaxis, intervalle, std)
      
      data_sp = data.frame(X=X ,Y=Y ,XerrorBar=XerrorBar ,YerrorBar=YerrorBar)
      table(intervalle)
      
      
      p9 = ggplot(data_sp,aes(x=X,y=Y,text=paste("Nb of samples by group",table(intervalle))))  + theme_bw() +
        ylab(names(which(axisIntra==input$y_intra))) + xlab(names(which(axisIntra==input$x_intra))) +
        geom_errorbar(aes(ymin=Y-YerrorBar, ymax=Y+YerrorBar),size=0.1) +
        geom_errorbarh(aes(xmin=X-XerrorBar, xmax=X+XerrorBar),size=0.1)+
        ggtitle(paste("No introns studied=",nrow(species_intron))) + 
        geom_point(pch=21,size=5,fill=color_species) + theme(
          axis.title.x = element_text(color="black", size=25,family="economica"),
          axis.title.y = element_text(color="black", size=25, family="economica"),
          axis.text.y =  element_text(color="black", size=20, family="economica"),
          axis.text.x =  element_text(color="black", size=20, family="economica"),
          title =  element_text(color="black", size=15, family="economica"),
          legend.text =  element_text(color="black", size=20, family="economica")
        )
    } else {
      
      data_sp = data.frame(X=unlist(listeAxis[input$x_intra]))
      
      p9 = ggplot(data_sp,aes(x=X)) + geom_histogram( position="dodge",bins=200,alpha=0.8,fill=color_species,col="black")+ theme_bw() +
        xlab(names(which(axisIntra==input$x_intra))) +ggtitle(paste("N=",nrow(data_sp))) + ylab("Density") + theme(
          axis.title.x = element_text(color="black", size=25,family="economica"),
          axis.title.y = element_text(color="black", size=25, family="economica"),
          axis.text.y =  element_text(color="black", size=20, family="economica"),
          axis.text.x =  element_text(color="black", size=20, family="economica"),
          title =  element_text(color="black", size=15, family="economica"),
          legend.text =  element_text(color="black", size=20, family="economica")
        )
    }
    
    if ( "Y log10" %in% input$scale_intra  ){ 
      p9 = p9 + scale_y_log10()
    }
    
    if ( "X log10" %in% input$scale_intra  ){
      p9 = p9 + scale_x_log10()
    } 
    
    ggplotly(p9,height = 700, width = 1100) 
  })
  
  ### ONGLET STRUCTURE GENE
  observe({
    if (input$tabs == "Gene structure"){
      species = input$species_gene_struct
      species = str_replace_all(species," ","_")
      species_genes = read.delim(paste("www/per_species_data/",species,"/by_gene_analysis.tab",sep=""), header=T , sep="\t",comment.char = "#")
      
      if ( grepl("busco_id_",input$gene_list) ){
        domain = str_replace(input$gene_list,"busco_id_","")
        
        busco_gene = read.delim(paste("www/per_species_data/",species,"/busco_to_gene_id_",domain,sep=""))
        rownames(busco_gene) = busco_gene$gene_id
        
        species_genes = species_genes[species_genes$gene_id %in% busco_gene$gene_id, ]
        species_genes$busco_id = busco_gene[species_genes$gene_id,]$busco_id
        nameGene <<- species_genes$gene_id
        names(nameGene) <<- paste(species_genes$busco_id,species_genes$gene_name,sep=" | ")
        
      } else {
        nameGene <<- species_genes$gene_id
        names(nameGene) <<- paste(species_genes$gene_id,species_genes$gene_name,sep=" | ")
      }
      updateSelectizeInput(session, "studied_gene", choices=nameGene,selected = nameGene[1],
                           options=list(maxOptions=100), server = T)
    }
  })
  
  output$species_image <- renderImage({
    secies_name <- str_replace(input$species_gene_struct," ","_")
    list(src=paste("www/species_images/",
                   secies_name,".png",sep=""))
  },deleteFile=FALSE)
  
  output$structureGene <- renderPlotly({
    species = input$species_gene_struct
    species = str_replace_all(species," ","_")
    id_selected = input$studied_gene
    if (id_selected != ""){
      if (id_selected != ""){
        start.time <- Sys.time()
        header = readLines(paste("www/per_species_data/",species,"/by_intron_cds.tab",sep=""),n=16)
        header = read.table(text=header)
        intron = system(paste("grep '",id_selected,"\t' ","www/per_species_data/",species,"/by_intron_cds.tab",sep=""),intern=T)
        print(intron)
        intron = read.table(text=intron)
        colnames(intron) = header
        intron$length = abs(intron$splice5 - intron$splice3)
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        print(time.taken)
        as.data.frame(intron,sep="\t")
        
        intron$category_intron = paste("Class:",intron$intron_class," in CDS:",intron$into_cds,sep="" )
        
        intron$position = paste("Sp3:",intron$splice3,"Sp5:",intron$splice5)
        intron$sum_n1 = as.numeric(intron$n1)
        p2 = ggplot(intron,aes(group=category_intron)) +
          geom_rect( aes(label=position , xmin=splice3,ymin=sum_n1*1.1,xmax=splice5,ymax=sum_n1,fill=category_intron),
                     size=0.5,alpha=1,col="black" ) +
          geom_segment(aes(x=splice3,y=sum_n1*1.1,xend=splice5+length/2*(splice3-splice5)/abs(splice5-splice3),
                           yend=sum_n1*1.4,fill=category_intron),
                       size=0.5,alpha=1,col="black")+
          geom_segment(aes(x=splice5,y=sum_n1*1.1,xend=splice3-length/2*(splice3-splice5)/abs(splice5-splice3),
                           yend=sum_n1*1.4,fill=category_intron),
                       size=0.5,alpha=1,col="black") +
          theme_bw()  +
          ylab("Sum of n1") +
          scale_fill_manual(name="Intron\ngroup",values=set_color) +
          xlab("Position on chromosome (bp)") +
          ggtitle(paste("Chromosome:",intron[1,"seqname"],"and Strand:",intron[1,"strand"])) +  theme(
            axis.title.x = element_text(color="black", size=25,family="economica"),
            axis.title.y = element_text(color="black", size=25, family="economica"),
            axis.text.y =  element_text(color="black", size=20, family="economica"),
            axis.text.x =  element_text(color="black", size=20, family="economica"),
            title =  element_text(color="black", size=15, family="economica"),
            text =  element_text(color="black", size=20, family="economica"),
            legend.text =  element_text(color="black", size=20, family="economica")
          ) + scale_y_log10()  +
          geom_vline(aes(xintercept=splice3,col=category_intron),alpha=0.1,fill="black") +
          geom_vline(aes(xintercept=splice5,col=category_intron),alpha=0.1,fill="black")+
          guides(col = guide_legend(override.aes = list( size = 6,alpha=1),
                                    label.theme = element_text(color="black",
                                                               size=26,face="italic", family="economica",vjust = 1.5,margin = margin(t = 5))))
        
        if (length(seq(min(intron$splice3,intron$splice5),
                       max(intron$splice3,intron$splice5),
                       input$sliderscale)) < 10){
          p2 = p2 + scale_x_continuous(breaks = seq(min(intron$splice3,intron$splice5),
                                                    max(intron$splice3,intron$splice5),
                                                    input$sliderscale) )
        }
        
        ggplotly(p2,height = 1000 * .8, )
      }
    }
  })
  
  ## PHYLOGENETIC TREE
  ranges2 <- reactiveValues(x = NULL, y = NULL, text_size=NULL)
  
  observe({
    brush <- input$plot2_brush
    if (!is.null(brush)) {
      ranges2$x <- c(brush$xmin, brush$xmax)
      ranges2$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges2$x <- NULL
      ranges2$y <- NULL
    }
  })
  
  observe({
    tree_name <- input$select_tree
    tree <- read.tree(tree_name)
    # tree$tip.label[grepl("dingo",tree$tip.label)] = "Canis_lupus"
    tree$tip.label <- str_replace_all(tree$tip.label,"_"," ")
    edge_group <- str_replace_all(tree$tip.label,"_"," ")
    edge_clade <- rep("branch",length(tree$edge[,2]))
    for (group in unique(edge_group)){
      if (group %in% unlist(listNomSpecies)){
        edge_clade[tree$edge[,2] %in% grep(group,edge_group)] =
          names(listNomSpecies[unlist(lapply(listNomSpecies,function(x) group %in% x))])
      }
    }
    
    for (clade in unique(edge_clade)){
      edge_clade[ which.edge(tree,  tree$edge[,2][edge_clade == clade] ) ] = clade
    }
    node_metadata = data.frame(node=tree$edge[,2],color=edge_clade)
    
    p = ggtree(tree, layout=input$layout_tree,size=1)  
    p <- p %<+% node_metadata  + aes(color=color) + 
      scale_color_manual("Clade",values=Clade_color[unique(edge_clade)]) +    theme(
        panel.background = element_rect(fill = "#f5f5f5", linetype = "dashed")
      ) 
    tree_plot <<- p
  })
  
  output$plot_principal <- renderPlot({
    input$layout_tree
    input$select_tree
    tree_plot + geom_tiplab(size=input$tip_size,nudge_x = input$spacing,colour="black")
  })
  
  output$plot_zoomed <- renderPlot({
    input$layout_tree
    input$select_tree
    tree_plot + geom_tiplab(size=input$tip_size_zoom,nudge_x = input$spacing_zoom,colour="black") +
      coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE)
  })
  
  
  output$download_fpkm <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), sep="")
    },
    content = function(file) {
      species = str_replace(input$species_selected_intra," ","_")
      data = read.delim(paste('www/per_species_data/',species,"/by_gene_analysis.tab",sep=""))
      write.table(data, file=file,row.names=F, col.names=T, sep="\t", quote=F)
    }
  )
  output$download_busco_id <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), sep="")
    },
    content = function(file) {
      species = str_replace(input$species_selected_intra," ","_")
      domain =input$busco_intra
      data = read.delim(paste("www/per_species_data/",species,"/busco_to_gene_id_",domain,sep="",sep=""))
      write.table(data, file=file,row.names=F, col.names=T, sep="\t", quote=F)
    }
  )
  
  output$download_svr <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), sep="")
    },
    content = function(file) {
      species = str_replace(input$species_selected_intra," ","_")
      data = read.delim(paste('www/per_species_data/',species,"/by_intron_major_overlap.tab",sep=""))
      write.table(data, file=file,row.names=F, col.names=T, sep="\t", quote=F)
    }
  )
}

# Complete app with UI and server components
shinyApp(ui, server)
