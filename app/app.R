library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(shinydisconnect)
library(shinycssloaders)
library(shinyjs)
library(reticulate)
library(DT)
library(readr)
library(shinyBS)
library(ggpubr)
library(shinyWidgets)
library(fenr)
library(stringr)
library(reticulate)

reticulate::install_python(version = '3.9')


virtualenv_create(envname = "fibrosis_shiny",
                  version  = "3.9")
conda_list()
virtualenv_list()



use_virtualenv("fibrosis_shiny")

py_install("scanpy")
py_install("python-igraph")
py_install("leidenalg")
py_install("decoupler")
py_install("omnipath")
py_install("marsilea")
py_install("pydeseq2")
py_install("adjustText")
py_install("psutil")

use_virtualenv("fibrosis_shiny")
sc <- import("scanpy")
dc <- import("decoupler")
pydeseq2_dds <- import("pydeseq2.dds")
pydeseq2_ds <- import("pydeseq2.ds")

ui <- fluidPage(
  disconnectMessage(
    text = "An error occurred. Please refresh the page and try again.",
    refresh = "Refresh",
    background = "#FFFFFF",
    colour = "#444444",
    refreshColour = "#2c3e50",
    overlayColour = "#000000",
    overlayOpacity = 0.6,
    width = 450,
    top = 50,
    size = 22,
    css = ""
  ),
  navbarPage(
    title = "Multi-species scRNA-seq Atlas of MASLD",
    theme = bs_theme(preset = 'flatly', base_font = 'Lato', code_font = 'Lato', heading_font = 'Lato'),
    tags$head(
      tags$link(rel = "icon", type = "image/png", sizes = "32x32", href = "tabicon.png")
    ),
    tabPanel(
      title = "Workflow",
      value = "tab_workflow",
      icon = icon("network-wired")
    ),
    tabPanel(
      title = "Explore & Analyze Datasets",
      value = "tab_explore_datasets",
      icon = icon("hourglass"),
      tabsetPanel(
        type = "pills",
        tabPanel(
          title = "Import Dataset",
          value = "subtab_import_dataset",
          icon = icon("file-import"),
          hr(),
          sidebarLayout(
            sidebarPanel(width = 2,
                         selectInput("selection_organism", "Select Organism",
                                     choices = c("Human", "Mouse", "Zebrafish", "Integrated")),
                         uiOutput("dataset_selection_list"),
                         actionButton("import_dataset", class = "btn-primary", "Load Dataset", width = '100%')
            ),
            mainPanel(
              fluidRow(
                column(width = 6, textOutput("textoutput_data_numbers")),
                column(width = 6, uiOutput("umap_col_selection"))
              ),
              fluidRow(
                column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_UMAP"))),
                column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_UMAP_selection")))
              ),
              br(),
              hr(),
              br(),
              fluidRow(
                column(width = 6, uiOutput("ranks_col_selection")),
                column(width = 6)
              ),
              fluidRow(
                column(width = 6, shinycssloaders::withSpinner(DTOutput("table_cluster_markers"))),
                column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_CellType_groups")))
              )
            )
          )
        ),
        tabPanel(
          title = "Cluster Selection",
          value = "subtab_cluster_selection",
          icon = icon("crosshairs"),
          hr(),
          sidebarLayout(
            sidebarPanel(width = 2,
                         bsTooltip(id = "selected_culusters",
                                   title = "This is your tooltip text",
                                   placement = "right",
                                   trigger = "hover"),
                         uiOutput("cluster_selection_list"),
                         uiOutput("button_cluster_filter"),
                         hr(),
                         selectInput("cluster_selection_visualization_type", "Select Visualization Method",
                                     choices = c("Visualize Expression of Gene",
                                                 "Visualize Expression of Geneset",
                                                 "Calculate Co-Expression")),
                         hr(),
                         conditionalPanel(condition = "input.cluster_selection_visualization_type == 'Visualize Expression of Gene'",
                                          uiOutput("gene_list_expression")),
                         conditionalPanel(condition = "input.cluster_selection_visualization_type == 'Visualize Expression of Geneset'",
                                          textInput("name_first_geneset", "Name First Geneset", placeholder = "Enter Name", width = '100%'),
                                          textAreaInput("first_geneset_list", "Enter Gene Names", placeholder = "Enter Names", width = '100%', height = '100%')

                                          ),
                         
                         uiOutput("gene_list_coexpresion_first"),
                         uiOutput("gene_list_coexpresion_second"),
                         uiOutput("button_visualize_cluster_selection"),
                         conditionalPanel(condition = "input.cluster_selection_visualization_type == 'Visualize Expression of Geneset'",
                                          br(),
                                          hr(),
                                          br(),
                                          uiOutput("geneset_name_input_second"),
                                          uiOutput("geneset_list_input_second"),
                                          uiOutput("geneset_list_go_second"))
            ),
            mainPanel(
              conditionalPanel(
                condition = "input.cluster_selection_visualization_type == 'Visualize Expression of Gene'",
                fluidRow(
                  column(width = 6, textOutput("info_cluster_selection")),
                  column(width = 6, uiOutput("umap_col_selection_clusters"))
                ),
                fluidRow(
                  column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_UMAP_expression_clusters"))),
                  column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_UMAP_selection_clusters")))
                ),
                br(),
                hr(),
                br(),
                fluidRow(
                  column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_violin_expression_clusters_groups"))),
                  column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_violin_expression_clusters")))
                )
              ),
              conditionalPanel(condition = "input.cluster_selection_visualization_type == 'Visualize Expression of Geneset'  && input.visualize_expression_genest != 0",
                               fluidRow(column(width = 6), column(width = 6,
                                                                  selectInput("first_geneset_enrichment_violin_selection", label = NULL, choices = c("Clusters", "Groups")))),
                               fluidRow(column(width = 6,
                                               shinycssloaders::withSpinner(imageOutput("image_output_enrichment_first_set"))), 
                                        column(width = 6,
                                               conditionalPanel(condition = "input.first_geneset_enrichment_violin_selection == 'Clusters'",
                                                                shinycssloaders::withSpinner(imageOutput("image_output_enrichment_first_violin_clusters"))
                                                                ),
                                               conditionalPanel(condition = "input.first_geneset_enrichment_violin_selection == 'Groups'",
                                                                shinycssloaders::withSpinner(imageOutput("image_output_enrichment_first_violin_groups"))
                                                                )
                                               )
                                               ),
                               br(),
                               hr(),
                               br(),
                               fluidRow(column(width = 6), column(width = 6,
                                                                  conditionalPanel(condition = "input.second_geneset_run != 0",
                                                                                   selectInput("second_geneset_enrichment_violin_selection", label = NULL, choices = c("Clusters", "Groups")))
                                                                  
                                                                  
                                                                  )),
                               fluidRow(column(width = 6,
                                               shinycssloaders::withSpinner(imageOutput("image_output_enrichment_second_set"))
                                               ), column(width = 6,
                                                         conditionalPanel(condition = "input.second_geneset_enrichment_violin_selection == 'Clusters'",
                                                                          shinycssloaders::withSpinner(imageOutput("image_output_enrichment_second_violin_clusters"))
                                                         ),
                                                         conditionalPanel(condition = "input.second_geneset_enrichment_violin_selection == 'Groups'",
                                                                          shinycssloaders::withSpinner(imageOutput("image_output_enrichment_second_violin_groups"))
                                                         )
                                                         
                                                         )))

              ,
              conditionalPanel(
                condition = "input.cluster_selection_visualization_type == 'Calculate Co-Expression' && input.visualize_cluster_selection != 0",
                fluidRow(
                  column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_UMAP_coexpression_first"))),
                  column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_UMAP_coexpression_second")))
                ),
                br(),
                hr(),
                br(),
                fluidRow(
                  column(width = 6,
                         radioGroupButtons(
                           inputId = "test_choice",
                           label = NULL,
                           choices = c("Spearman", "Pearson"),
                           justified = TRUE, width = '100%'
                         ),
                         bslib::input_switch("remove_zero_counts", "Remove Zero Counts", width = '100%')
                  ),
                  column(width = 6, uiOutput("top_correlated_genes"))
                ),
                fluidRow(
                  column(width = 6, shinycssloaders::withSpinner(plotOutput("stats_cor_plot"))),
                  column(width = 6,
                         fluidRow(
                           column(width = 6, shinycssloaders::withSpinner(DTOutput("first_gene_correlation_table"))),
                           column(width = 6, shinycssloaders::withSpinner(DTOutput("second_gene_correlation_table")))
                         )
                  )
                )
              )
            )
          )
        ),
        tabPanel(
          title = "Differential Expression",
          value = "subtab_differential_expression",
          icon = icon("bars-staggered"),
          hr(),
          sidebarLayout(
            sidebarPanel(width = 2,
                         radioGroupButtons(
                           inputId = "de_data",
                           label = NULL,
                           choices = c("All Data", "Filtered Data"),
                           justified = TRUE, width = '100%'
                         ),
                         radioGroupButtons(
                           inputId = "de_type",
                           label = NULL,
                           choices = c("Clusters", "Groups"),
                           justified = TRUE, width = '100%'
                         ),
                         hr(),
                         uiOutput("de_ident_selection_name_first"),
                         uiOutput("de_ident_selection_first"),
                         uiOutput("de_ident_selection_name_second"),
                         uiOutput("de_ident_selection_second"),
                         uiOutput("de_method_selection"),
                         uiOutput("de_run_analysis")
            ),
            mainPanel(
              conditionalPanel(
                condition = "input.de_run_dge > 0",
                fluidRow(
                  column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_dge_ranks"))),
                  column(width = 6, uiOutput("dge_table_filter"),
                         shinycssloaders::withSpinner(DTOutput("dge_dt"))
                  )
                ),
                br(),
                hr(),
                br(),
                fluidRow(
                  column(width = 6, shinycssloaders::withSpinner(imageOutput("imageoutput_dge_violin"))),
                  column(width = 6,
                         radioGroupButtons(
                           inputId = "de_radiobutton",
                           label = NULL,
                           choices = c("Violin Plot", "Enrichment"),
                           justified = TRUE, width = '100%'
                         ),
                         conditionalPanel(
                           condition = "input.de_radiobutton == 'Violin Plot'",
                           shinycssloaders::withSpinner(imageOutput("violin_dge"))
                         ),
                         conditionalPanel(
                           condition = "input.de_radiobutton == 'Enrichment'",
                           selectInput("de_enrichment_type", label = NULL, choices = c("GO", "BP", "KEGG", "Reactome", "WikiPathways")),
                           shinycssloaders::withSpinner(DTOutput("de_enrichment_table"))
                         )
                  )
                )
              )
            )
          )
        ),
        tabPanel(
          title = "Pseudo Bulk",
          value = "subtab_pseudo_buk",
          icon = icon("bucket"),
          hr(),
          sidebarLayout(
            sidebarPanel(width = 2,
                         radioGroupButtons(inputId = "pseudo_bulk_data",
                                           label = NULL,
                                           choices = c("All Data", "Filtered Data"),
                                           justified = TRUE, width = '100%'
                         ),
                         uiOutput("pseudo_ident_selection_first"),
                         uiOutput("pseudo_ident_selection_second"),
                         uiOutput("run_pseudo_bulk"),
                         hr(),
                         uiOutput("pdata_clusters"),
                         uiOutput("run_deseq2_ui")
            ),
            mainPanel(
              conditionalPanel(
                condition = "input.run_pseudo_bulk_act > 0",
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("pseudo_pca_selection_ui")
                  ),
                  column(width = 6,
                         uiOutput("selection_pseudo_volcano_table_ui"))
                ),
                fluidRow(
                  column(width = 6, shinycssloaders::withSpinner(imageOutput("pca_pseudo_bulk"))),
                  column(width = 6,
                         conditionalPanel(condition = "input.selection_pseudo_volcano_table == 'Volcano'",
                                          shinycssloaders::withSpinner(imageOutput("pca_pseudo_bulk_volcano"))),
                         conditionalPanel(condition = "input.selection_pseudo_volcano_table == 'Result Table'",
                                          shinycssloaders::withSpinner(DTOutput("pca_pseudo_bulk_results_table"))))
                ),
                br(),
                hr(),
                br(),
                fluidRow(
                  column(width = 6, uiOutput("pseudo_ora_selection_ui"),
                         
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'CollecTRI'",
                                          uiOutput("collectri_outputs")),
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'PROGENy'",
                                          uiOutput("progeny_outputs")),
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'MSigDB Hallmark'",
                                          uiOutput("msigdb_outputs"))
                         ),
                  column(width = 6, uiOutput("enrichment_1_pseudo"))
                ),
                fluidRow(
                  column(width = 6,
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'CollecTRI'",
                                          conditionalPanel(condition = "input.collectri_plot_selection == 'Top 25'",
                                                           shinycssloaders::withSpinner(imageOutput("collectri_top"))),
                                          conditionalPanel(condition = "input.collectri_plot_selection == 'Volcano'",
                                                           shinycssloaders::withSpinner(imageOutput("collectri_volcano"))),
                                          conditionalPanel(condition = "input.collectri_plot_selection == 'Network'",
                                                           div(
                                                             style = "overflow-x: auto; overflow-y: auto;",
                                                             shinycssloaders::withSpinner(imageOutput("collectri_network"))
                                                           ))),
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'PROGENy'",
                                          conditionalPanel(condition = "input.progeny_plot_selection == 'Scores'",
                                                           shinycssloaders::withSpinner(imageOutput("progeny_top"))),
                                          conditionalPanel(condition = "input.progeny_plot_selection == 'Targets'",
                                                           shinycssloaders::withSpinner(imageOutput("progeny_targets")))
                                          
                                          ),
                         conditionalPanel(condition = "input.pseudo_ora_selection == 'MSigDB Hallmark'",
                                          conditionalPanel(condition = "input.msigdb_plot_selection == 'Top Signatures'",
                                                           shinycssloaders::withSpinner(imageOutput("msigdb_top"))),
                                          conditionalPanel(condition = "input.msigdb_plot_selection == 'Running Score'",
                                                           shinycssloaders::withSpinner(imageOutput("msigdb_running")))
                                          
                         )
                         ),
                  column(width = 6, shinycssloaders::withSpinner(DTOutput("pseudo_enrichment_table")))
                )
              )
            )
          )
        )
        # ,
        # tabPanel(
        #   title = "Cell-Cell Interaction Analysis",
        #   value = "subtab_cell_cell_interaction_analysis_1",
        #   icon = icon("phone-volume")
        # )
      )
    )
    # ,
    # tabPanel(
    #   title = "Datasets",
    #   value = "tab_datasets",
    #   icon = icon("database")
    # ),
    # tabPanel(
    #   title = "About",
    #   value = "tab_about",
    #   icon = icon("info")
    # )
  )
)




server <- function(input, output,session) {
  
  ##################################### Import Dataset
  
  
  # Import Dataset
  ##ui actions
  # Read the dataset mapping CSV
  dataset_mapping <- read.csv("./dataset_mapping.csv", stringsAsFactors = FALSE)
  output$dataset_selection_list <- renderUI({
  # Filter datasets based on the selected organism
  filtered_datasets <- dataset_mapping %>%
    filter(species == input$selection_organism) %>%
    pull(dataset)
  
  # Check if there are any datasets for the selected species
  if (length(filtered_datasets) > 0) {
    selectInput("selection_dataset", "Select Dataset", choices = filtered_datasets)
  } else {
    selectInput("selection_dataset", "Select Dataset", choices = c("No datasets available"))
  }
})

  
  
  
  adata <- eventReactive(input$import_dataset, {
    # Import
    adata <- sc$read_h5ad(paste0("datasets/", input$selection_organism, "/", input$selection_dataset, ".h5ad"))
    # Disable the button
    shinyjs::disable("import_dataset")
    updateActionButton(session, "import_dataset", label = "", icon("circle-check"))
    # Return the loaded dataset
    return(adata)
  })
  
  
  #ui outputs
    output$textoutput_data_numbers <- renderText({
      req(adata())
      paste0("Number of Cells: ", adata()$n_obs,"   Number of Genes: ", adata()$n_vars)
    })
  
    output$imageoutput_UMAP <- renderImage({
    req(adata())
    sc$pl$umap(adata(),color = c('CellType'), legend_loc = "on data", show=FALSE,save = '.png')
    list(src = "figures/umap.png")
    }, deleteFile = TRUE)
    
    
    
    
    output$umap_col_selection <- renderUI({
      req(adata())
      # Use the column names as choices in selectInput
      selectInput("selection_umap_identity_select", "Select Identity", choices = colnames(adata()$obs))
    })
    
    output$imageoutput_UMAP_selection <- renderImage({
      req(adata(),input$selection_umap_identity_select)
      sc$pl$umap(adata(),color = input$selection_umap_identity_select, show=FALSE, save = 'second.png')
      list(src = "figures/umapsecond.png")
    }, deleteFile = TRUE)
  
    
    
    output$ranks_col_selection <- renderUI({
      req(adata())
      # Use the column names as choices in selectInput
      cluster_names <- unique(adata()$obs$CellType)
      cluster_names <- sort(cluster_names)
      cluster_names <- as.character(cluster_names)
      selectInput("selection_rank_select", "Select Cluster", choices = cluster_names)
    })
    
    
    

    output$imageoutput_CellType_groups <- renderImage({
      req(adata(),input$selection_rank_select)
      group_name <- list(input$selection_rank_select)
      # Convert the group number to a string
      sc$pl$rank_genes_groups(adata(), groups = group_name, sharey = FALSE, show=FALSE, save = '.png')
      list(src = "figures/rank_genes_groups_CellType.png")
    }, deleteFile = TRUE)
    
    

    
    output$table_cluster_markers <- renderDT({
      req(adata(),input$selection_rank_select)
      data_table_dges <- sc$get$rank_genes_groups_df(adata(), group= input$selection_rank_select)
      datatable(
        data_table_dges, 
        options = list(
          pageLength = 5,
          autoWidth = TRUE,
          scrollX = TRUE,
          scrollY = 300,  # Set the height for the scroll
          scroller = TRUE,
          dom = 'l<"toolbar">rtip',  # Removed 'f' from dom
          deferRender = TRUE
        ),
        filter = 'top',
        rownames = FALSE,
        extensions = 'Scroller'
      )%>%
        formatRound(which(sapply(data_table_dges, is.numeric)), digits = 3)
    })
    
    
    ##################################### Cluster Selection
    
    output$cluster_selection_list <- renderUI({
      req(adata())
      
      choices_for_selector <- unique(adata()$obs$CellType)
      choices_for_selector <- sort(choices_for_selector)
      choices_for_selector <- as.character(choices_for_selector) 
      selectizeInput(
        "selected_culusters", 
          "Choose Cluster",
        choices = choices_for_selector, 
        multiple = TRUE, options = list(placeholder = "Select Clusters to Filter"))
    })
    
    
    output$button_cluster_filter <- renderUI({
      req(adata(),input$selected_culusters)
      actionButton("filter_dataset_cluster_selection", "Filter Dataset", class = "btn-primary", icon = icon("filter"), width = '100%')
      
    })

    
    filtered_adata <- eventReactive(input$filter_dataset_cluster_selection, {

      selected_culusters <- list(input$selected_culusters)
      # filtered_adata <- adata()[adata()$obs['CellType'] == selected_culusters]
      
      filtered_adata = adata()[adata()$obs$CellType %in% input$selected_culusters, ]
      
      return(filtered_adata)
    })
    
    
    gene_list_adata <- eventReactive(c(input$import_dataset,input$filter_dataset_cluster_selection), {
      req(adata())
      
      if(is.null(input$filter_dataset_cluster_selection)){
        expression_table <- adata()$var
        gene_list_adata <- rownames(expression_table)
      }else{
        expression_table <- filtered_adata()$var
        gene_list_adata <- rownames(expression_table)
      }
      return(gene_list_adata)
    })
    
    output$gene_list_expression <- renderUI({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(gene_list_adata())
        gene_list_adata <- sort(gene_list_adata())
      selectizeInput(
        "gene_selection_cluster_expression", "Select Gene", choices = gene_list_adata
      )}
    })
    
    output$button_visualize_cluster_selection <- renderUI({
      req(gene_list_adata())
      actionButton("visualize_cluster_selection", "Visualize",icon = icon("eye"), width = '100%', class = "btn-primary")
    })
    
    
    output$umap_col_selection_clusters <- renderUI({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$visualize_cluster_selection)
      # Use the column names as choices in selectInput
      selectInput("selection_umap_identity_select_clusters", "Select Identity", choices = colnames(adata()$obs),
                  selected = 'CellType')}
    })
    
    
    output$info_cluster_selection <- renderText({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$visualize_cluster_selection)
      
      if(is.null(input$filter_dataset_cluster_selection)){
        paste0("Number of Cells: ", adata()$n_obs,"   Number of Genes: ", adata()$n_vars)
      }else{
        paste0("Number of Cells: ", filtered_adata()$n_obs,"   Number of Genes: ", filtered_adata()$n_vars)
      }
        }
      
      
    })
    
    
    output$imageoutput_UMAP_selection_clusters <- renderImage({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$selection_umap_identity_select_clusters,input$visualize_cluster_selection,input$gene_selection_cluster_expression)
      if(is.null(input$filter_dataset_cluster_selection)){
        sc$pl$umap(adata(),color = input$selection_umap_identity_select_clusters, show=FALSE, save = 'clusters.png')
        list(src = "figures/umapclusters.png")
      }else{
        sc$pl$umap(filtered_adata(),color = input$selection_umap_identity_select_clusters, show=FALSE, save = 'clusters.png')
        list(src = "figures/umapclusters.png")
      }
      }
      
    }, deleteFile = TRUE)
    
    output$imageoutput_UMAP_expression_clusters <- renderImage({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$visualize_cluster_selection,input$gene_selection_cluster_expression)
      if(is.null(input$filter_dataset_cluster_selection)){
        sc$pl$umap(adata(),color = input$gene_selection_cluster_expression, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'clusters_exp.png')
        list(src = "figures/umapclusters_exp.png")
      }else{
        sc$pl$umap(filtered_adata(),color = input$gene_selection_cluster_expression, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'clusters_exp.png')
        list(src = "figures/umapclusters_exp.png")
      }
      }
    }, deleteFile = TRUE)
    
    output$imageoutput_violin_expression_clusters <- renderImage({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$visualize_cluster_selection,input$gene_selection_cluster_expression)
      if(is.null(input$filter_dataset_cluster_selection)){
        sc$pl$violin(adata(), keys = input$gene_selection_cluster_expression, groupby = 'CellType', use_raw=F, layer = 'scvi_normalized', show=FALSE, rotation=90, save = "violin_exp.png")
        list(src = "figures/violinviolin_exp.png")
      }else{
        sc$pl$violin(filtered_adata(), keys = input$gene_selection_cluster_expression, groupby = 'CellType', use_raw=F, layer = 'scvi_normalized', show=FALSE, rotation=90, save = 'violin_exp.png')
        list(src = "figures/violinviolin_exp.png")
      }
      }
    }, deleteFile = TRUE)
    
    
    output$imageoutput_violin_expression_clusters_groups <- renderImage({
      if(input$cluster_selection_visualization_type == 'Visualize Expression of Gene'){
        
      req(input$visualize_cluster_selection,input$gene_selection_cluster_expression)
      if(is.null(input$filter_dataset_cluster_selection)){
        sc$pl$violin(adata(), keys = input$gene_selection_cluster_expression, groupby = 'Group', use_raw=F, layer = 'scvi_normalized', show=FALSE, rotation=90, save = "clusters_violin_exp.png")
        list(src = "figures/violinclusters_violin_exp.png")
      }else{
        sc$pl$violin(filtered_adata(), keys = input$gene_selection_cluster_expression, groupby = 'Group', use_raw=F, layer = 'scvi_normalized', show=FALSE, rotation=90, save = 'clusters_violin_exp.png')
        list(src = "figures/violinclusters_violin_exp.png")
      }
      }
    }, deleteFile = TRUE)
    
    ################################################################ Visualize genesets
    

    observeEvent(input$first_geneset_list, {
      genes <- input$first_geneset_list
      
      # Remove special characters except for "-"
      genes <- str_replace_all(genes, "[^A-Za-z0-9\\-\\s,]", "")
      
      # Split by various separators
      gene_list <- str_split(genes, "\\s+|,|\\t")[[1]]
      
      # Remove empty strings and duplicates, and trim whitespace
      gene_list <- unique(trimws(gene_list[gene_list != ""]))
      
      if(is.null(input$filter_dataset_cluster_selection)){
        expression_table <- adata()$raw$var
        gene_list_adata <- rownames(expression_table)
      }else{
        expression_table <- filtered_adata()$raw$var
        gene_list_adata <- rownames(expression_table)
      }
      
      gene_list <- intersect(gene_list, gene_list_adata)
      
      # Create a cleaned string with one gene name per line
      cleaned_genes <- paste(gene_list, collapse = "\n")
      

      # Update the text area input with the cleaned gene names
      
      updateTextAreaInput(session, "first_geneset_list", value = cleaned_genes)
    })
    
    
    first_gene_set_calc <- eventReactive(input$visualize_cluster_selection,{
      req(input$cluster_selection_visualization_type== "Visualize Expression of Geneset", input$first_geneset_list)
      
      gene_list <- strsplit(input$first_geneset_list, split = "\n")[[1]]
      
      if(!is.null(input$name_first_geneset)){
        net <- data.frame(
          source = rep(input$name_first_geneset, length(gene_list)),
          target = c(gene_list)
        )
      }else{
        net <- data.frame(
          source = rep("First Geneset", length(gene_list)),
          target = c(gene_list)
        )
      }

      if(is.null(input$filter_dataset_cluster_selection)){
        dc$run_ulm(adata(), net, weight=NULL, min_n=0)
        acts = dc$get_acts(adata(), 'ulm_estimate')
        acts_matrix <- acts$X
        acts_v <- as.vector(acts$X)
        max_e <- max(acts_v[is.finite(acts_v)], na.rm = TRUE)
        acts_matrix[!is.finite(acts_matrix)] <- max_e
        acts$X <- acts_matrix
        acts
        
        return(acts)
        
      }else{
        dc$run_ulm(filtered_adata(), net, weight=NULL, min_n=0)
        acts = dc$get_acts(filtered_adata(), 'ulm_estimate')
        acts_matrix <- acts$X
        acts_v <- as.vector(acts$X)
        max_e <- max(acts_v[is.finite(acts_v)], na.rm = TRUE)
        acts_matrix[!is.finite(acts_matrix)] <- max_e
        acts$X <- acts_matrix
        acts
  
        return(acts)
        
      }
    
    })
    
    
    output$image_output_enrichment_first_set <- renderImage({
      req(input$visualize_cluster_selection,input$first_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_first_geneset)){
        sc$pl$umap(first_gene_set_calc(), color= rownames(first_gene_set_calc()$var), cmap='RdBu_r', show = F, vmax = 5, save = "first_gene_enrichment_umap.png")
      }else{
        sc$pl$umap(first_gene_set_calc(), color= rownames(first_gene_set_calc()$var), cmap='RdBu_r', show = F, vmax = 5, save = "first_gene_enrichment_umap.png")
      }
      list(src = "figures/umapfirst_gene_enrichment_umap.png")
      
    }, deleteFile = TRUE)
    
    
    output$image_output_enrichment_first_violin_clusters <- renderImage({
      req(input$visualize_cluster_selection,input$first_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_first_geneset)){
        sc$pl$violin(first_gene_set_calc(), keys = rownames(first_gene_set_calc()$var), groupby='CellType', show = F, rotation=90, save = "first_gene_enrichment_violin.png")
      }else{
        sc$pl$violin(first_gene_set_calc(), keys = rownames(first_gene_set_calc()$var), groupby='CellType', show = F, rotation=90, save = "first_gene_enrichment_violin.png")
      }
      list(src = "figures/violinfirst_gene_enrichment_violin.png")
      
    }, deleteFile = TRUE)
    
    output$image_output_enrichment_first_violin_groups <- renderImage({
      req(input$visualize_cluster_selection,input$first_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_first_geneset)){
        sc$pl$violin(first_gene_set_calc(), keys = rownames(first_gene_set_calc()$var), groupby='Group', show = F, rotation=90, save = "first_gene_enrichment_violin_group.png")
      }else{
        sc$pl$violin(first_gene_set_calc(), keys = rownames(first_gene_set_calc()$var), groupby='Group', show = F, rotation=90, save = "first_gene_enrichment_violin_group.png")
      }
      list(src = "figures/violinfirst_gene_enrichment_violin_group.png")
      
    }, deleteFile = TRUE)
    
    
    output$geneset_name_input_second <- renderUI({
      req(input$visualize_cluster_selection)
      textInput("name_second_geneset", "Name second Geneset", placeholder = "Enter Name", width = '100%')
    })
    
    output$geneset_list_input_second <- renderUI({
      req(input$visualize_cluster_selection)
      textAreaInput("second_geneset_list", "Enter Gene Names", placeholder = "Enter Names", width = '100%')
    })
    
    output$geneset_list_go_second <- renderUI({
      req(input$second_geneset_list)
      actionButton("second_geneset_run", "Visualize", width = '100%', icon = icon("eye"), class = "btn-primary")
    })
    
    observeEvent(input$second_geneset_list, {
      genes <- input$second_geneset_list
      
      # Remove special characters except for "-"
      genes <- str_replace_all(genes, "[^A-Za-z0-9\\-\\s,]", "")
      
      # Split by various separators
      gene_list <- str_split(genes, "\\s+|,|\\t")[[1]]
      
      # Remove empty strings and duplicates, and trim whitespace
      gene_list <- unique(trimws(gene_list[gene_list != ""]))
      
      gene_list <- intersect(gene_list, gene_list_adata())
      
      # Create a cleaned string with one gene name per line
      cleaned_genes <- paste(gene_list, collapse = "\n")
      
      
      # Update the text area input with the cleaned gene names
      
      updateTextAreaInput(session, "second_geneset_list", value = cleaned_genes)
    })
    
    
    
    
    
    second_gene_set_calc <- eventReactive(input$visualize_cluster_selection,{
      req(input$cluster_selection_visualization_type== "Visualize Expression of Geneset", input$second_geneset_run)
      
      gene_list <- strsplit(input$second_geneset_list, split = "\n")[[1]]
      
      if(!is.null(input$name_second_geneset)){
        net <- data.frame(
          source = rep(input$name_second_geneset, length(gene_list)),
          target = gene_list
        )
      }else{
        net <- data.frame(
          source = rep("Second Geneset", length(gene_list)),
          target = gene_list
        )
      }
      
      
      
      if(is.null(input$filter_dataset_cluster_selection)){
        dc$run_ulm(adata(), net, weight=NULL, min_n=0)
        acts = dc$get_acts(adata(), 'ulm_estimate')
        acts_matrix <- acts$X
        acts_v <- as.vector(acts$X)
        max_e <- max(acts_v[is.finite(acts_v)], na.rm = TRUE)
        acts_matrix[!is.finite(acts_matrix)] <- max_e
        acts$X <- acts_matrix
        acts
        
        return(acts)
        
      }else{
        dc$run_ulm(filtered_adata(), net, weight=NULL, min_n=0)
        acts = dc$get_acts(filtered_adata(), 'ulm_estimate')
        acts_matrix <- acts$X
        acts_v <- as.vector(acts$X)
        max_e <- max(acts_v[is.finite(acts_v)], na.rm = TRUE)
        acts_matrix[!is.finite(acts_matrix)] <- max_e
        acts$X <- acts_matrix
        acts
        
        return(acts)
        
      }
      
    })
    
    
    output$image_output_enrichment_second_set <- renderImage({
      req(input$second_geneset_run,input$second_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_second_geneset)){
        sc$pl$umap(second_gene_set_calc(), color= rownames(second_gene_set_calc()$var), cmap='RdBu_r', show = F,vmax = 5, save = "second_gene_enrichment_umap.png")
      }else{
        sc$pl$umap(second_gene_set_calc(), color= rownames(second_gene_set_calc()$var), cmap='RdBu_r', show = F,vmax = 5, save = "second_gene_enrichment_umap.png")
      }
      list(src = "figures/umapsecond_gene_enrichment_umap.png")
      
    }, deleteFile = TRUE)
    
    
    output$image_output_enrichment_second_violin_clusters <- renderImage({
      req(input$second_geneset_run,input$second_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_second_geneset)){
        sc$pl$violin(second_gene_set_calc(), keys = rownames(second_gene_set_calc()$var), groupby='CellType', show = F, rotation=90, save = "second_gene_enrichment_violin.png")
      }else{
        sc$pl$violin(second_gene_set_calc(), keys = rownames(second_gene_set_calc()$var), groupby='CellType', show = F, rotation=90, save = "second_gene_enrichment_violin.png")
      }
      list(src = "figures/violinsecond_gene_enrichment_violin.png")
      
    }, deleteFile = TRUE)
    
    output$image_output_enrichment_second_violin_groups <- renderImage({
      req(input$second_geneset_run,input$second_geneset_list, input$cluster_selection_visualization_type == 'Visualize Expression of Geneset')
      if(!is.null(input$name_second_geneset)){
        sc$pl$violin(second_gene_set_calc(), keys = rownames(second_gene_set_calc()$var), groupby='Group', show = F, rotation=90, save = "second_gene_enrichment_violin_group.png")
      }else{
        sc$pl$violin(second_gene_set_calc(), keys = rownames(second_gene_set_calc()$var), groupby='Group', show = F, rotation=90, save = "second_gene_enrichment_violin_group.png")
      }
      list(src = "figures/violinsecond_gene_enrichment_violin_group.png")
      
    }, deleteFile = TRUE)
    
    
    
    
    
    
    
    
    
    
    ## co-expression
    
    output$gene_list_coexpresion_first <- renderUI({
      if(input$cluster_selection_visualization_type == 'Calculate Co-Expression'){
        
        req(gene_list_adata())
        gene_list_adata <- sort(gene_list_adata())
        
        selectizeInput(
          "gene_selection_cluster_coexpression_first", "Select First Gene", choices = gene_list_adata
        )}
    })
    
    output$gene_list_coexpresion_second <- renderUI({
      if(input$cluster_selection_visualization_type == 'Calculate Co-Expression'){
        
        req(gene_list_adata())
        gene_list_adata <- sort(gene_list_adata())
        selectizeInput(
          "gene_selection_cluster_coexpression_second", "Select Second Gene", choices = gene_list_adata
        )}
    })
    
    output$imageoutput_UMAP_coexpression_first <- renderImage({
      if(input$cluster_selection_visualization_type == 'Calculate Co-Expression'){
        
        req(input$visualize_cluster_selection, input$gene_selection_cluster_coexpression_first)
        if(is.null(input$filter_dataset_cluster_selection)){
          sc$pl$umap(adata(),color = input$gene_selection_cluster_coexpression_first, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'coexp_1.png')
          list(src = "figures/umapcoexp_1.png")
        }else{
          sc$pl$umap(filtered_adata(),color = input$gene_selection_cluster_coexpression_first, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'coexp_1.png')
          list(src = "figures/umapcoexp_1.png")
        }
      }
    }, deleteFile = TRUE)
    
    output$imageoutput_UMAP_coexpression_second <- renderImage({
      if(input$cluster_selection_visualization_type == 'Calculate Co-Expression'){
        
        req(input$visualize_cluster_selection,input$gene_selection_cluster_coexpression_second)
        if(is.null(input$filter_dataset_cluster_selection)){
          sc$pl$umap(adata(),color = input$gene_selection_cluster_coexpression_second, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'coexp_2.png')
          list(src = "figures/umapcoexp_2.png")
        }else{
          sc$pl$umap(filtered_adata(),color = input$gene_selection_cluster_coexpression_second, layer = 'scvi_normalized', vmax = 5, show=FALSE, save = 'coexp_2.png')
          list(src = "figures/umapcoexp_2.png")
        }
      }
    }, deleteFile = TRUE)
    
    statistics_coexpression <- reactive({
      req(input$gene_selection_cluster_coexpression_first,input$gene_selection_cluster_coexpression_second)
      if(is.null(input$filter_dataset_cluster_selection)){
        counts <- adata()$X
        first_gene_count <- counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_first)]
        second_gene_count <- counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_second)]
      }else{
        counts <- filtered_adata()$X
        first_gene_count <- counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_first)]
        second_gene_count <- counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_second)]
      }
      
      
      
      if(input$gene_selection_cluster_coexpression_first == input$gene_selection_cluster_coexpression_second){
        first_gene <- input$gene_selection_cluster_coexpression_first
        second_gene = paste0(input$gene_selection_cluster_coexpression_second,"_")
      }else{
        first_gene = input$gene_selection_cluster_coexpression_first
        second_gene = input$gene_selection_cluster_coexpression_second
      }
      statistics_coexpression <- cbind(first_gene_count, second_gene_count)
      statistics_coexpression <- as.data.frame(statistics_coexpression)
      colnames(statistics_coexpression) <- c(first_gene, second_gene)
      results <- list(statistics_coexpression,first_gene, second_gene)
      return(results)
    })
    
    # output$spearman_results <- renderText({
    #   req(statistics_coexpression(), input$visualize_cluster_selection, input$gene_selection_cluster_coexpression_first,input$gene_selection_cluster_coexpression_second)
    #   paste(paste0(input$gene_selection_cluster_coexpression_first, "-", input$gene_selection_cluster_coexpression_second),
    #     paste0("Test Statistic: ", statistics_coexpression()[[1]]$statistic),
    #         paste0("Spearman's rho: ", statistics_coexpression()[[1]]$estimate),
    #         paste0("p-value: ", statistics_coexpression()[[1]]$p.value),
    #         collapse = "\n")
    #   
    # })
    
    output$stats_cor_plot <- renderPlot({
      req(statistics_coexpression(), input$visualize_cluster_selection)
      if(input$remove_zero_counts == F & input$test_choice == "Spearman"){
        ggscatter(data = statistics_coexpression()[[1]],
                  x = statistics_coexpression()[[2]], 
                  y = statistics_coexpression()[[3]],
                  title = "Spearman Correlation",
                  add = "reg.line", 
                  conf.int = TRUE,                add.params = list(color = "2c3e50",
                                                                    fill = "lightgray"))+
          stat_cor(method = "spearman", label.x = 0.5, label.sep = "\n", size = 10, color = "red")
      }else if(input$remove_zero_counts == T && input$test_choice == "Spearman"){
        data_for_plot <- statistics_coexpression()[[1]]
        data_for_plot <- data_for_plot %>%
          filter_all(all_vars(. != 0))
        
        ggscatter(data = data_for_plot,
                  x = statistics_coexpression()[[2]], 
                  y = statistics_coexpression()[[3]],
                  title = "Spearman Correlation",
                  add = "reg.line", 
                  conf.int = TRUE,
                  add.params = list(color = "2c3e50", fill = "lightgray")) +
          stat_cor(method = "spearman", label.x = 0.5, label.sep = "\n", size = 10, color = "red")
        
      }else if(input$remove_zero_counts == F && input$test_choice == "Pearson"){
        ggscatter(data = statistics_coexpression()[[1]],
                  x = statistics_coexpression()[[2]], 
                  y = statistics_coexpression()[[3]],
                  title = "Pearson Correlation",
                  add = "reg.line", 
                  conf.int = TRUE,                add.params = list(color = "2c3e50",
                                                                    fill = "lightgray"))+
          stat_cor(method = "pearson", label.x = 0.5, label.sep = "\n", size = 10, color = "red")
      }else if(input$remove_zero_counts == T && input$test_choice == "Pearson"){
        data_for_plot <- statistics_coexpression()[[1]]
        data_for_plot <- data_for_plot %>%
          filter_all(all_vars(. != 0))
        ggscatter(data = data_for_plot,
                  x = statistics_coexpression()[[2]], 
                  y = statistics_coexpression()[[3]],
                  title = "Pearson Correlation",
                  add = "reg.line", 
                  conf.int = TRUE,
                  add.params = list(color = "2c3e50", fill = "lightgray")) +
          stat_cor(method = "pearson", label.x = 0.5, label.sep = "\n", size = 10, color = "red")
      }
    })
    
    ## calculate top correlated genes
    
    
    output$top_correlated_genes <- renderUI({
      req(input$visualize_cluster_selection, statistics_coexpression())
      fluidRow(column(width = 6,
                      actionButton("top_correlated_first_gene",paste0("Find Correlated Genes with ", input$gene_selection_cluster_coexpression_first), icon = icon("chart-line"),class = "btn-primary", width = '100%')),
               column(width = 6,
                      actionButton("top_correlated_second_gene",paste0("Find Correlated Genes with ", input$gene_selection_cluster_coexpression_second), icon = icon("chart-line"),class = "btn-primary", width = '100%')))
      
    })
    
    
    correlation_table_first_gene <- eventReactive(input$top_correlated_first_gene,{
      if(is.null(input$filter_dataset_cluster_selection)){
        
      normalized_counts <- as.matrix(adata()$X)
      normalized_counts <- as.data.frame(normalized_counts)
      colnames(normalized_counts) <- gene_list_adata()
      first_gene_count <- normalized_counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_first)]
      
      correlation_df <- sapply(names(normalized_counts), function(x) {
        test_result <- cor.test(first_gene_count, normalized_counts[[x]], method = "spearman")
        return(c(correlation = test_result$estimate, p_value = test_result$p.value))
      })
      
      correlation_df <- as.data.frame(t(correlation_df))
      num_tests <- ncol(normalized_counts)
      correlation_df$Bonferroni_p_value <- pmin(1, correlation_df$p_value * num_tests)
      correlation_df$Gene <- rownames(correlation_df)
      colnames(correlation_df)[1:2] <- c("Correlation","p-val")
      correlation_df <- correlation_df[, c("Gene", "Correlation", "p-val", "Bonferroni_p_value")]
      
      return(correlation_df)}else{
        normalized_counts <- as.matrix(filtered_adata()$X)
        normalized_counts <- as.data.frame(normalized_counts)
        colnames(normalized_counts) <- gene_list_adata()
        first_gene_count <- normalized_counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_first)]
        
        correlation_df <- sapply(names(normalized_counts), function(x) {
          test_result <- cor.test(first_gene_count, normalized_counts[[x]], method = "spearman")
          return(c(correlation = test_result$estimate, p_value = test_result$p.value))
        })
        
        correlation_df <- as.data.frame(t(correlation_df))
        num_tests <- ncol(normalized_counts)
        correlation_df$Bonferroni_p_value <- pmin(1, correlation_df$p_value * num_tests)
        correlation_df$Gene <- rownames(correlation_df)
        colnames(correlation_df)[1:2] <- c("Correlation","p-val")
        correlation_df <- correlation_df[, c("Gene", "Correlation", "p-val", "Bonferroni_p_value")]
        
        return(correlation_df)
        
      }
    })
    

    output$first_gene_correlation_table <- renderDT({
      req(correlation_table_first_gene())
      datatable(
        correlation_table_first_gene(), 
        options = list(
          pageLength = 5,
          autoWidth = TRUE,
          scrollX = TRUE),
        filter = 'top',
        rownames = FALSE,
        extensions = 'Scroller'
      )%>%
        formatRound(columns=c("Correlation","p-val","Bonferroni_p_value"), digits=3)
      
    })
    
    
    correlation_table_second_gene <- eventReactive(input$top_correlated_second_gene,{
      if(is.null(input$filter_dataset_cluster_selection)){
        
      normalized_counts <- as.data.frame(adata()$X)
      colnames(normalized_counts) <- gene_list_adata()
      second_gene_count <- normalized_counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_second)]
      
      correlation_df <- sapply(names(normalized_counts), function(x) {
        test_result <- cor.test(second_gene_count, normalized_counts[[x]], method = "spearman")
        return(c(correlation = test_result$estimate, p_value = test_result$p.value))
      })
      
      correlation_df <- as.data.frame(t(correlation_df))
      num_tests <- ncol(normalized_counts)
      correlation_df$Bonferroni_p_value <- pmin(1, correlation_df$p_value * num_tests)
      correlation_df$Gene <- rownames(correlation_df)
      colnames(correlation_df)[1:2] <- c("Correlation","p-val")
      correlation_df <- correlation_df[, c("Gene", "Correlation", "p-val", "Bonferroni_p_value")]
      
      return(correlation_df)}else{
        normalized_counts <- as.data.frame(filtered_adata()$X)
        colnames(normalized_counts) <- gene_list_adata()
        second_gene_count <- normalized_counts[,which(gene_list_adata() == input$gene_selection_cluster_coexpression_second)]
        
        correlation_df <- sapply(names(normalized_counts), function(x) {
          test_result <- cor.test(second_gene_count, normalized_counts[[x]], method = "spearman")
          return(c(correlation = test_result$estimate, p_value = test_result$p.value))
        })
        
        correlation_df <- as.data.frame(t(correlation_df))
        num_tests <- ncol(normalized_counts)
        correlation_df$Bonferroni_p_value <- pmin(1, correlation_df$p_value * num_tests)
        correlation_df$Gene <- rownames(correlation_df)
        colnames(correlation_df)[1:2] <- c("Correlation","p-val")
        correlation_df <- correlation_df[, c("Gene", "Correlation", "p-val", "Bonferroni_p_value")]
        
        
        return(correlation_df)
      }
    })
    
    
    output$second_gene_correlation_table <- renderDT({
      req(correlation_table_second_gene())
      datatable(
        correlation_table_second_gene(), 
        options = list(
          pageLength = 5,
          autoWidth = FALSE,
          scrollX = TRUE),
        filter = 'top',
        rownames = FALSE,
        extensions = 'Scroller'
      ) %>%
        formatRound(columns=c("Correlation","p-val","Bonferroni_p_value"), digits=3)
      
    })
    
    ##########################################################
    
    ######### Differential Expression
    
    output$de_ident_selection_name_first <- renderUI({
      req(adata())
      if(input$de_data == "All Data" && input$de_type == "Clusters"){
        textInput("de_ident_1_name", "Name First Ident", placeholder = "Enter Name")
      }else if(input$de_data == "All Data" && input$de_type == "Groups"){
        textInput("de_ident_1_name", "Name First Ident", placeholder = "Enter Name")
      }else if(input$de_data == "Filtered Data" && input$de_type == "Clusters"){
        if(!is.null(filtered_adata())){
          textInput("de_ident_1_name", "Name First Ident", placeholder = "Enter Name")
        }else if(is.null(filtered_adata())){
          "You Should First Filter Data"
        }
      }else if(input$de_data == "Filtered Data" && input$de_type == "Groups"){
        if(!is.null(filtered_adata())){
          textInput("de_ident_1_name", "Name First Ident", placeholder = "Enter Name")
        }else{
          "You Should First Filter Data"
        }
      }
    })
    
    
    output$de_ident_selection_first <- renderUI({
      req(adata(),input$de_ident_1_name)
      if(input$de_data == "All Data" && input$de_type == "Clusters"){
        selectizeInput("de_ident_1", "First Ident", choices = sort(unique(adata()$obs$CellType)), multiple = T)
      }else if(input$de_data == "All Data" && input$de_type == "Groups"){
        selectizeInput("de_ident_1", "First Ident", choices = sort(unique(adata()$obs$Group)), multiple = T)
      }else if(input$de_data == "Filtered Data" && input$de_type == "Clusters"){
        if(!is.null(filtered_adata())){
          selectizeInput("de_ident_1", "First Ident", choices = sort(unique(filtered_adata()$obs$CellType)), multiple = T)
        }else if(is.null(filtered_adata())){
          "You Should First Filter Data"
        }
      }else if(input$de_data == "Filtered Data" && input$de_type == "Groups"){
        if(!is.null(filtered_adata())){
          selectizeInput("de_ident_1", "First Ident", choices = sort(unique(filtered_adata()$obs$Group)), multiple = T)
        }else{
          "You Should First Filter Data"
        }
      }
    })
    
    
    output$de_ident_selection_name_second <- renderUI({
      req(adata(), input$de_ident_1)
      if(input$de_data == "All Data" && input$de_type == "Clusters"){
        textInput("de_ident_2_name", "Name Second Ident", placeholder = "Enter Name")
      }else if(input$de_data == "All Data" && input$de_type == "Groups"){
        textInput("de_ident_2_name", "Name Second Ident", placeholder = "Enter Name")
      }else if(input$de_data == "Filtered Data" && input$de_type == "Clusters"){
        if(!is.null(filtered_adata())){
          textInput("de_ident_2_name", "Name Second Ident", placeholder = "Enter Name")
        }else if(is.null(filtered_adata())){
          "You Should First Filter Data"
        }
      }else if(input$de_data == "Filtered Data" && input$de_type == "Groups"){
        if(!is.null(filtered_adata())){
          textInput("de_ident_2_name", "Name Second Ident", placeholder = "Enter Name")
        }else{
          "You Should First Filter Data"
        }
      }
    })
    
    
    
    output$de_ident_selection_second <- renderUI({
      req(adata(),input$de_ident_2_name)
      if(input$de_data == "All Data" && input$de_type == "Clusters"){
        first_ident <- input$de_ident_1
        idents <- unique(adata()$obs$CellType)
        first_ident <- idents %in% first_ident
        idents <- idents[!first_ident]
        selectizeInput("de_ident_2", "Second Ident", choices = sort(idents),multiple = T)
      }else if(input$de_data == "All Data" && input$de_type == "Groups"){
        first_ident <- input$de_ident_1
        idents <- unique(adata()$obs$Group)
        first_ident <- idents %in% first_ident
        idents <- idents[!first_ident]
        selectizeInput("de_ident_2", "Second Ident", choices = sort(idents),multiple = T)
      }else if(input$de_data == "Filtered Data" && input$de_type == "Clusters"){
        first_ident <- input$de_ident_1
        idents <- unique(filtered_adata()$obs$CellType)
        first_ident <- idents %in% first_ident
        idents <- idents[!first_ident]
        selectizeInput("de_ident_2", "Second Ident", choices = sort(idents),multiple = T)
      }else if(input$de_data == "Filtered Data" && input$de_type == "Groups"){
        first_ident <- input$de_ident_1
        idents <- unique(filtered_adata()$obs$Group)
        first_ident <- idents %in% first_ident
        idents <- idents[!first_ident]
        selectizeInput("de_ident_2", "Second Ident", choices = sort(idents),multiple = T)
      }
    })
    
    output$de_method_selection <- renderUI({
      req(input$de_ident_2)
      selectInput("de_method_selection", "Select Method", choices = c('t-test', 'wilcoxon', 't-test_overestim_var'), selected = "wilcoxon")
    })
    
    output$de_run_analysis <- renderUI({
      req(input$de_method_selection)
      actionButton("de_run_dge", "Run DGE", class = "btn-primary", icon = icon("play"), width = '100%')
    })
    
    
    de_filtering_identifying_clusters <- eventReactive(input$de_run_dge,{
      
      if(input$de_data == "All Data" && input$de_type == "Clusters"){
        adata <- adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$de_ident_1
        second_identity <- input$de_ident_2
        first_identity_rows <- obs_data$CellType %in% first_identity
        second_identity_rows <- obs_data$CellType %in% second_identity
        obs_data$ident[first_identity_rows] <- input$de_ident_1_name 
        obs_data$ident[second_identity_rows] <- input$de_ident_2_name       
        adata$obs['ident'] = obs_data$ident
        selection_idents <- list(input$de_ident_1_name, input$de_ident_2_name)
        adata = adata[adata$obs$ident %in% selection_idents, ]

        return(adata)
        

      }else if(input$de_data == "All Data" && input$de_type == "Groups"){
        adata <- adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$de_ident_1
        second_identity <- input$de_ident_2
        first_identity_rows <- obs_data$Group %in% first_identity
        second_identity_rows <- obs_data$Group %in% second_identity
        obs_data$ident[first_identity_rows] <- input$de_ident_1_name 
        obs_data$ident[second_identity_rows] <- input$de_ident_2_name       
        adata$obs$ident = obs_data$ident
        selection_idents <- c(input$de_ident_1_name, input$de_ident_2_name)
        adata = adata[adata$obs$ident %in% selection_idents, ]
        return(adata)
        
      }else if(input$de_data == "Filtered Data" && input$de_type == "Clusters"){
        adata <- filtered_adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$de_ident_1
        second_identity <- input$de_ident_2
        first_identity_rows <- obs_data$CellType %in% first_identity
        second_identity_rows <- obs_data$CellType %in% second_identity
        obs_data$ident[first_identity_rows] <- input$de_ident_1_name 
        obs_data$ident[second_identity_rows] <- input$de_ident_2_name       
        adata$obs['ident'] = obs_data$ident
        selection_idents <- c(input$de_ident_1_name, input$de_ident_2_name)
        adata = adata[adata$obs$ident %in% selection_idents, ]
        return(adata)
        
      }else if(input$de_data == "Filtered Data" && input$de_type == "Groups"){
        adata <- filtered_adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$de_ident_1
        second_identity <- input$de_ident_2
        first_identity_rows <- obs_data$Group %in% first_identity
        second_identity_rows <- obs_data$Group %in% second_identity
        obs_data$ident[first_identity_rows] <- input$de_ident_1_name 
        obs_data$ident[second_identity_rows] <- input$de_ident_2_name       
        adata$obs['ident'] = obs_data$ident
        selection_idents <- c(input$de_ident_1_name, input$de_ident_2_name)
        adata = adata[adata$obs$ident %in% selection_idents, ]
        
        return(adata)
      }
    })
    
    
    de_dge_calculation <- reactive({
      req(de_filtering_identifying_clusters())
        adata <- de_filtering_identifying_clusters()
        sc$tl$rank_genes_groups(adata, 'ident', groups=list(input$de_ident_1_name), reference= input$de_ident_2_name, method= input$de_method_selection,pts = T)
        return(adata)

    })
    
    
    output$imageoutput_dge_ranks <- renderImage({
      req(de_dge_calculation())
      group_name <- list(input$de_ident_1_name)
      sc$pl$rank_genes_groups(de_dge_calculation(), groups = group_name, sharey = FALSE, show=FALSE, save = 'dge.png')
      list(src = "figures/rank_genes_groups_identdge.png")
    }, deleteFile = TRUE)
  
    output$dge_dt <- renderDT({
      req(de_dge_calculation())
      group_name <- list(input$de_ident_1_name)
      result_df = sc$get$rank_genes_groups_df(de_dge_calculation(), group = group_name, key = 'rank_genes_groups')
      
      colnames(result_df) <- c("Gene", "Scores", "LogFC", "p-val", "adj-p", "pct")
      
      if(input$dge_filter_table_logfc > 0){
        result_df <- result_df %>%
          filter(LogFC > 2)
      }
      
      if(input$dge_filter_table_logfc_neg > 0){
        result_df <- result_df %>%
          filter(LogFC < -2)
      }
      
      if(input$dge_filter_table_score_neg > 0){
        result_df <- result_df %>%
          filter(Scores < -2)
      }
      
      if(input$dge_filter_table_score > 0){
        result_df <- result_df %>%
          filter(Scores > 2)
      }
      
      if(input$dge_filter_table_pval > 0){
        result_df <- result_df %>%
          filter(`adj-p` < 0.05)
      }
      
      if(input$dge_filter_table_reset > 0){
          result_df = sc$get$rank_genes_groups_df(de_dge_calculation(), group = group_name, key = 'rank_genes_groups')
          
          colnames(result_df) <- c("Gene", "Scores", "LogFC", "p-val", "adj-p", "pct")
          
      }
      
      datatable(
        result_df, 
        options = list(
          pageLength = 5,
          autoWidth = TRUE,
          scrollX = TRUE,
          dom = 'l<"toolbar">rtip' 
        ),
        filter = 'top',
        rownames = FALSE,
        selection = "single",
        extensions = 'Scroller'
      ) %>%
        formatRound(which(sapply(result_df, is.numeric)), digits = 3)
    })
    
    output$dge_table_filter <- renderUI({
      req(de_dge_calculation())
      fluidRow(
        column(width = 2, actionButton("dge_filter_table_logfc", label = "LogFC > 2", class = "btn-primary", width = '100%')),
        column(width = 2, actionButton("dge_filter_table_logfc_neg", label = "LogFC < -2", class = "btn-primary", width = '100%')),
        column(width = 2, actionButton("dge_filter_table_pval", label = "Adjusted-P Val < 0.05", class = "btn-primary", width = '100%')),
        column(width = 2, actionButton("dge_filter_table_reset", label = "Reset Filters", class = "btn-primary", width = '100%')),
        column(width = 2, actionButton("dge_filter_table_score", label = "Score > 2", class = "btn-primary", width = '100%')),
        column(width = 2, actionButton("dge_filter_table_score_neg", label = "Score < -2", class = "btn-primary", width = '100%'))
      )
    })
    
    
    
    
    output$imageoutput_dge_violin <- renderImage({
      req(de_dge_calculation())
      group_name <- list(input$de_ident_1_name)
      sc$pl$rank_genes_groups_violin(de_dge_calculation(), groups = group_name,show=FALSE, save = '.png')

      list(src = paste0("figures/rank_genes_groups_ident_",group_name,".png"))
      
    }, deleteFile = TRUE)
    
    
    output$radiobutton_dge <- renderUI({
      req(de_dge_calculation())
      fluidRow(column(width = 6, actionButton("dge_violin","Violin Plot", icon = icon("chart-simple"),class = "btn-primary", width = '100%')),
               column(width = 6, actionButton("dge_go","Go Enrichment", icon = icon("chart-bar"),class = "btn-primary", width = '100%')))
    })
    
    output$violin_dge <- renderImage({
      req(de_dge_calculation(), input$dge_dt_rows_selected)

      group_name <- list(input$de_ident_1_name)
      result_df = sc$get$rank_genes_groups_df(de_dge_calculation(), group = group_name, key = 'rank_genes_groups')
      selected_row_index <- input$dge_dt_rows_selected
      selected_row <- result_df[selected_row_index, ]
      selected_value <- selected_row[[1]]
      
      sc$pl$violin(de_dge_calculation(), keys = selected_value, groupby = 'ident', show=FALSE, rotation=90, save = "violin_dge.png")
      list(src = "figures/violinviolin_dge.png")
    }, deleteFile = TRUE)
    
    
    de_enrichment_calc <- reactive({
      req(de_dge_calculation(), input$dge_dt_rows_all)
      group_name <- list(input$de_ident_1_name)
      result_df = sc$get$rank_genes_groups_df(de_dge_calculation(), group = group_name, key = 'rank_genes_groups')

      if(input$selection_organism %in% c("Human", "Integrated")){
        load("enrichment_sets/human.RData")
        all_genes <- de_dge_calculation()$var
        all_genes <- rownames(all_genes)
        selected_genes <- result_df[input$dge_dt_rows_all, ]
        selected_genes <- selected_genes[,1]
        
        
        bp_human <- prepare_for_enrichment(bp_human$terms, bp_human$mapping, all_genes, feature_name = "gene_symbol")
        go_human <- prepare_for_enrichment(go_human$terms, go_human$mapping, all_genes, feature_name = "gene_symbol")
        kegg_human <- prepare_for_enrichment(kegg_human$terms, kegg_human$mapping, all_genes, feature_name = "gene_symbol")
        reactome_human <- prepare_for_enrichment(reactome_human$terms, reactome_human$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_human <- prepare_for_enrichment(wiki_pathways_human$terms, wiki_pathways_human$mapping, all_genes, feature_name = "text_label")
        
        bp_human <- functional_enrichment(all_genes, selected_genes, bp_human)
        go_human <- functional_enrichment(all_genes, selected_genes, go_human)
        kegg_human <- functional_enrichment(all_genes, selected_genes, kegg_human)
        reactome_human <- functional_enrichment(all_genes, selected_genes, reactome_human)
        wiki_pathways_human <- functional_enrichment(all_genes, selected_genes, wiki_pathways_human)
        
        de_enrichment_calc <- list(bp_human, go_human, kegg_human, reactome_human, wiki_pathways_human)
        
        return(de_enrichment_calc)
        
      }else if(input$selection_organism == "Mouse"){
        
        load("enrichment_sets/mouse.RData")
        all_genes <- de_dge_calculation()$var
        all_genes <- rownames(all_genes)
        selected_genes <- result_df[input$dge_dt_rows_all, ]
        selected_genes <- selected_genes[,1]
        
        bp_mouse <- prepare_for_enrichment(bp_mouse$terms, bp_mouse$mapping, all_genes, feature_name = "gene_symbol")
        go_mouse <- prepare_for_enrichment(go_mouse$terms, go_mouse$mapping, all_genes, feature_name = "gene_symbol")
        kegg_mouse <- prepare_for_enrichment(kegg_mouse$terms, kegg_mouse$mapping, all_genes, feature_name = "gene_symbol")
        reactome_mouse <- prepare_for_enrichment(reactome_mouse$terms, reactome_mouse$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_mouse <- prepare_for_enrichment(wiki_pathways_mouse$terms, wiki_pathways_mouse$mapping, all_genes, feature_name = "text_label")
        
        bp_mouse <- functional_enrichment(all_genes, selected_genes, bp_mouse)
        go_mouse <- functional_enrichment(all_genes, selected_genes, go_mouse)
        kegg_mouse <- functional_enrichment(all_genes, selected_genes, kegg_mouse)
        reactome_mouse <- functional_enrichment(all_genes, selected_genes, reactome_mouse)
        wiki_pathways_mouse <- functional_enrichment(all_genes, selected_genes, wiki_pathways_mouse)
        de_enrichment_calc <- list(bp_mouse, go_mouse, kegg_mouse, reactome_mouse, wiki_pathways_mouse)
        return(de_enrichment_calc)
      }else if(input$selection_organism == "Zebrafish"){
        
        load("enrichment_sets/zebrafish.RData")
        all_genes <- de_dge_calculation()$var
        all_genes <- rownames(all_genes)
        selected_genes <- result_df[input$dge_dt_rows_all, ]
        selected_genes <- selected_genes[,1]
        
        bp_zebrafish <- prepare_for_enrichment(bp_zebrafish$terms, bp_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        go_zebrafish <- prepare_for_enrichment(go_zebrafish$terms, go_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        kegg_zebrafish <- prepare_for_enrichment(kegg_zebrafish$terms, kegg_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        reactome_zebrafish <- prepare_for_enrichment(reactome_zebrafish$terms, reactome_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_zebrafish <- prepare_for_enrichment(wiki_pathways_zebrafish$terms, wiki_pathways_zebrafish$mapping, all_genes, feature_name = "text_label")
        
        bp_zebrafish <- functional_enrichment(all_genes, selected_genes, bp_zebrafish)
        go_zebrafish <- functional_enrichment(all_genes, selected_genes, go_zebrafish)
        kegg_zebrafish <- functional_enrichment(all_genes, selected_genes, kegg_zebrafish)
        reactome_zebrafish <- functional_enrichment(all_genes, selected_genes, reactome_zebrafish)
        wiki_pathways_zebrafish <- functional_enrichment(all_genes, selected_genes, wiki_pathways_zebrafish)
        de_enrichment_calc <- list(bp_zebrafish, go_zebrafish, kegg_zebrafish, reactome_zebrafish, wiki_pathways_zebrafish)
        return(de_enrichment_calc)
        
      }
      # else if(input$selection_organism == "Integrated"){
      #   load("enrichment_sets/human.RData")
      #   all_genes <- de_dge_calculation()$var
      #   all_genes <- rownames(all_genes)
      #   selected_genes <- result_df[input$dge_dt_rows_all, ]
      #   
      #   bp_human <- prepare_for_enrichment(bp_human$terms, bp_human$mapping, all_genes, feature_name = "gene_symbol")
      #   go_human <- prepare_for_enrichment(go_human$terms, go_human$mapping, all_genes, feature_name = "gene_symbol")
      #   kegg_human <- prepare_for_enrichment(kegg_human$terms, kegg_human$mapping, all_genes, feature_name = "gene_symbol")
      #   reactome_human <- prepare_for_enrichment(reactome_human$terms, reactome_human$mapping, all_genes, feature_name = "gene_symbol")
      #   wiki_pathways_human <- prepare_for_enrichment(wiki_pathways_human$terms, wiki_pathways_human$mapping, all_genes, feature_name = "gene_symbol")
      #   
      #   bp_human <- functional_enrichment(all_genes, selected_genes, bp_human)
      #   go_human <- functional_enrichment(all_genes, selected_genes, go_human)
      #   kegg_human <- functional_enrichment(all_genes, selected_genes, kegg_human)
      #   reactome_human <- functional_enrichment(all_genes, selected_genes, reactome_human)
      #   wiki_pathways_human <- functional_enrichment(all_genes, selected_genes, wiki_pathways_human)
      #   de_enrichment_calc <- list(bp_human, go_human, kegg_human, reactome_human, wiki_pathways_human)
      #   return(de_enrichment_calc)}
    })

    output$de_enrichment_table <- renderDT({
      req(de_enrichment_calc())
      if(input$de_enrichment_type == "GO"){
        datatable(
          de_enrichment_calc()[[2]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip' 
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(de_enrichment_calc()[[2]], is.numeric)), digits = 3)
      }else if(input$de_enrichment_type == "BP"){
        datatable(
          de_enrichment_calc()[[1]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip' 
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(de_enrichment_calc()[[1]], is.numeric)), digits = 3)
      }else if(input$de_enrichment_type == "KEGG"){
        datatable(
          de_enrichment_calc()[[3]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip' 
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(de_enrichment_calc()[[3]], is.numeric)), digits = 3)
      }else if(input$de_enrichment_type == "Reactome"){
        datatable(
          de_enrichment_calc()[[4]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip' 
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(de_enrichment_calc()[[4]], is.numeric)), digits = 3)
      }else if(input$de_enrichment_type == "WikiPathways"){
        datatable(
          de_enrichment_calc()[[5]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip' 
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(de_enrichment_calc()[[5]], is.numeric)), digits = 3)
      }
    })
    
  
    # selectizeInput("de_ident_2", "Second Ident", choices = sort(unique(adata()$obs$CellType)),multiple = T)
    
    
    ########## pseudo_bulk
    
    output$pseudo_ident_selection_first <- renderUI({
      req(adata())
      if(input$pseudo_bulk_data == "All Data"){
        selectizeInput("pseudo_ident_1", "First Ident", choices = sort(unique(adata()$obs$Group)), multiple = T)
      }else if(input$pseudo_bulk_data == "Filtered Data"){
        selectizeInput("pseudo_ident_1", "First Ident", choices = sort(unique(filtered_adata()$obs$Group)), multiple = T)
      }
    })
    
    output$pseudo_ident_selection_second <- renderUI({
      req(adata(),input$pseudo_ident_1)
      if(input$pseudo_bulk_data == "All Data"){
        idents <- unique(adata()$obs$Group)
        first_ident <- idents %in% input$pseudo_ident_1
        idents <- idents[!first_ident]
        selectizeInput("pseudo_ident_2", "Second Ident", choices = sort(unique(idents)), multiple = T)
      }else if(input$pseudo_bulk_data == "Filtered Data"){
        idents <- unique(filtered_adata()$obs$Group)
        first_ident <- idents %in% input$pseudo_ident_1
        idents <- idents[!first_ident]
        selectizeInput("pseudo_ident_2", "Second Ident", choices = sort(unique(selectizeInput)), multiple = T)
      }
    })
    
    output$run_pseudo_bulk <- renderUI({
      req(input$pseudo_ident_2)
      actionButton("run_pseudo_bulk_act", "Run Pseudo Bulk", class = "btn-primary",width = '100%')
    })
    
    output$pseudo_pca_selection_ui  <- renderUI({
      req(input$run_pseudo_bulk_act)
      fluidRow(column(width = 6, selectInput("pseudo_pca_selection", label = NULL, choices = c("ident", "CellType"))),
               column(width = 6, actionButton("associations_pseudo_modal", "Show PC Associations", class = "btn-primary", width = '100%')))
    })
    
    
    
    pdata <- eventReactive(input$run_pseudo_bulk_act, {
      if(input$pseudo_bulk_data == "All Data"){
        adata <- adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$pseudo_ident_1
        second_identity <- input$pseudo_ident_2
        first_identity_rows <- obs_data$Group %in% first_identity
        second_identity_rows <- obs_data$Group %in% second_identity
        obs_data$ident[first_identity_rows] <- 'Selection' 
        obs_data$ident[second_identity_rows] <- 'Reference'       
        adata$obs['ident'] = obs_data$ident
        selection_idents <- list('Selection', 'Reference' )
        adata = adata[adata$obs$ident %in% selection_idents, ]  
      }else if(input$pseudo_bulk_data == "Filtered Data"){
        adata <- filtered_adata()
        obs_data <- adata$obs
        obs_data$ident <- "ident"
        first_identity <- input$pseudo_ident_1
        second_identity <- input$pseudo_ident_2
        first_identity_rows <- obs_data$Group %in% first_identity
        second_identity_rows <- obs_data$Group %in% second_identity
        obs_data$ident[first_identity_rows] <- 'Selection' 
        obs_data$ident[second_identity_rows] <- 'Reference'       
        adata$obs['ident'] = obs_data$ident
        selection_idents <- list('Selection', 'Reference' )
        adata = adata[adata$obs$ident %in% selection_idents, ]
      }
      
      pdata <- dc$get_pseudobulk(
        adata,
        sample_col = 'Group',
        groups_col = 'CellType',
        layer = 'counts',
        mode = 'sum',
        min_cells = as.integer(10),
        min_counts = as.integer(1000)
      )
      pdata$layers['counts'] <- pdata$X
      sc$pp$normalize_total(pdata, target_sum = 1e4)
      sc$pp$log1p(pdata)
      sc$pp$scale(pdata, max_value = 10)
      sc$tl$pca(pdata)
      dc$swap_layer(pdata, 'counts', X_layer_key = "None", inplace = TRUE)
      dc$get_metadata_associations(
        pdata,
        obs_keys = list('ident', 'CellType', 'psbulk_n_cells', 'psbulk_counts'),
        obsm_key = 'X_pca',
        uns_key = 'pca_anova',
        inplace = TRUE
      )
      return(pdata)
    })
    
    output$pca_pseudo_bulk <- renderImage({
      req(pdata())
      if(input$pseudo_pca_selection == "ident"){
        sc$pl$pca(pdata(), color = c('ident'), show = F, save = "pseudo_pca.png")
        list(src = "figures/pcapseudo_pca.png")
      }else if(input$pseudo_pca_selection == "CellType"){
        sc$pl$pca(pdata(), color = c('CellType'), show = F, save = "pseudo_pca.png")
        list(src = "figures/pcapseudo_pca.png")
      }else if(is.null(input$pseudo_pca_selection)){
        sc$pl$pca(pdata(), color = c('ident'), show = F, save = "pseudo_pca.png")
        list(src = "figures/pcapseudo_pca.png")
      }
    }, deleteFile = TRUE)
    
    output$pca_pseudo_bulk_associations <- renderImage({
      req(pdata())
      dc$plot_associations(
        pdata(),
        uns_key = 'pca_anova',
        obsm_key = 'X_pca',
        stat_col = 'p_adj',
        obs_annotation_cols = list('ident', 'CellType'),
        titles = c('Principle component scores', 'Adjusted p-values from ANOVA'),
        figsize = c(as.integer(7), as.integer(5)),
        n_factors = as.integer(10),
        save = "adjusted_pca.png"
      )
      list(src = "adjusted_pca.png")
    }, deleteFile = TRUE)
    
    
    observeEvent(input$associations_pseudo_modal, {
      showModal(modalDialog(
        div(
          style = "height: 500px; overflow-y: auto;",
          imageOutput("pca_pseudo_bulk_associations")
        ),
        size = "l", 
        easyClose = TRUE
      ))
    })
    
    output$pdata_clusters <- renderUI({
      req(pdata())
      selectizeInput("pdata_clusters_filter", "Select Clusters to Filter", choices = sort(unique(pdata()$obs$CellType)),
                     multiple = T)
    })
    
    output$run_deseq2_ui <- renderUI({
      req(input$pdata_clusters_filter)
      actionButton("run_deseq2", "Run DESeq2", class = "btn-primary",width = '100%')
    })
    
    results_df <- eventReactive(input$run_deseq2,{
      
      filtered_pseudo <- pdata()[pdata()$obs$CellType %in% input$pdata_clusters_filter, ]
      genes <- dc$filter_by_expr(filtered_pseudo, group = 'ident', min_count = 10, min_total_count = 15)
      filtered_pseudo <- filtered_pseudo[, genes]
      filtered_pseudo <- filtered_pseudo$copy()
      
      dds <- pydeseq2_dds$DeseqDataSet(
        adata = filtered_pseudo,
        design_factors = 'ident',
        ref_level=list('ident', 'Reference'),
        refit_cooks = TRUE
      )
      
      dds$deseq2()
      
      stat_res = pydeseq2_ds$DeseqStats(
        dds,
        contrast=list("ident", 'Selection', 'Reference')
      )
      
      
      stat_res$summary()
      
      results_df = stat_res$results_df
      results_df <- as.data.frame(results_df)
      results_df[["Gene Name"]] <- rownames(results_df)
      results_df <- results_df %>% relocate(`Gene Name`, .before = 1)
      mat <- as.data.frame(t(results_df['stat']))
      rownames(mat) <- "Selected_Clusters"
      results_df <- list(results_df, mat)
      return(results_df)
      
    })
    
    output$selection_pseudo_volcano_table_ui <- renderUI({
      req(results_df()[[1]][[1]])
      selectInput("selection_pseudo_volcano_table", label=NULL, choices = c("Result Table", "Volcano"))
    })
    
    
    output$pca_pseudo_bulk_volcano <- renderImage({
      req(results_df()[[1]][[1]])
      
      dataset <- results_df()[[1]]
      dataset[is.na(dataset)] <- 0.000001
      dataset[dataset == ""] <- 0.000001
      
      dc$plot_volcano_df(
        results_df()[[1]],
        x='log2FoldChange',
        y='padj',
        top= as.integer(20),
        save = "pseudo_volcano.png"
      )
      list(src = "pseudo_volcano.png")
    }, deleteFile = TRUE)
    
    output$pca_pseudo_bulk_results_table <- renderDT({
      req(results_df()[[1]])
      dataset <- results_df()[[1]]
      dataset[is.na(dataset)] <- 0
      dataset[dataset == ""] <- 0
      
      datatable(
        dataset, 
        options = list(
          pageLength = 5,
          autoWidth = TRUE,
          scrollX = TRUE,
          dom = 'l<"toolbar">rtip'  # Removed 'f' from dom
        ),
        filter = 'top',
        rownames = FALSE,
        extensions = 'Scroller'
      )
    })
    
    output$enrichment_1_pseudo <- renderUI({
      req(results_df())
      fluidRow(column(width = 6,
                      selectInput("pseudo_enrichment_type", label = NULL, choices = c("GO", "BP", "KEGG", "Reactome", "WikiPathways"))),
               column(width = 6,
                      actionButton("run_pseudo_enrichment", class = "btn-primary", "Calculate", width = '100%')
                      ))
      
      
    })
    
    
    pseudo_enrichment_calc <- eventReactive(input$run_pseudo_enrichment,{
      req(results_df(), input$pca_pseudo_bulk_results_table_rows_all, input$pseudo_enrichment_type)
      
      if(input$selection_organism %in% c("Human", "Integrated")){
        load("enrichment_sets/human.RData")
        all_genes <- rownames(results_df()[[1]])
        
        selected_genes <- rownames(results_df()[[1]][input$pca_pseudo_bulk_results_table_rows_all, ])
        
        
        
        bp_human <- prepare_for_enrichment(bp_human$terms, bp_human$mapping, all_genes, feature_name = "gene_symbol")
        go_human <- prepare_for_enrichment(go_human$terms, go_human$mapping, all_genes, feature_name = "gene_symbol")
        kegg_human <- prepare_for_enrichment(kegg_human$terms, kegg_human$mapping, all_genes, feature_name = "gene_symbol")
        reactome_human <- prepare_for_enrichment(reactome_human$terms, reactome_human$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_human <- prepare_for_enrichment(wiki_pathways_human$terms, wiki_pathways_human$mapping, all_genes, feature_name = "text_label")
        
        bp_human <- functional_enrichment(all_genes, selected_genes, bp_human)
        go_human <- functional_enrichment(all_genes, selected_genes, go_human)
        kegg_human <- functional_enrichment(all_genes, selected_genes, kegg_human)
        reactome_human <- functional_enrichment(all_genes, selected_genes, reactome_human)
        wiki_pathways_human <- functional_enrichment(all_genes, selected_genes, wiki_pathways_human)
        
        de_enrichment_calc <- list(bp_human, go_human, kegg_human, reactome_human, wiki_pathways_human)
        
        return(de_enrichment_calc)
        
      }else if(input$selection_organism == "Mouse"){
        
        load("enrichment_sets/mouse.RData")
        all_genes <- rownames(results_df()[[1]])
        
        selected_genes <- rownames(results_df()[[1]][input$pca_pseudo_bulk_results_table_rows_all, ])
        
        bp_mouse <- prepare_for_enrichment(bp_mouse$terms, bp_mouse$mapping, all_genes, feature_name = "gene_symbol")
        go_mouse <- prepare_for_enrichment(go_mouse$terms, go_mouse$mapping, all_genes, feature_name = "gene_symbol")
        kegg_mouse <- prepare_for_enrichment(kegg_mouse$terms, kegg_mouse$mapping, all_genes, feature_name = "gene_symbol")
        reactome_mouse <- prepare_for_enrichment(reactome_mouse$terms, reactome_mouse$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_mouse <- prepare_for_enrichment(wiki_pathways_mouse$terms, wiki_pathways_mouse$mapping, all_genes, feature_name = "text_label")
        
        bp_mouse <- functional_enrichment(all_genes, selected_genes, bp_mouse)
        go_mouse <- functional_enrichment(all_genes, selected_genes, go_mouse)
        kegg_mouse <- functional_enrichment(all_genes, selected_genes, kegg_mouse)
        reactome_mouse <- functional_enrichment(all_genes, selected_genes, reactome_mouse)
        wiki_pathways_mouse <- functional_enrichment(all_genes, selected_genes, wiki_pathways_mouse)
        de_enrichment_calc <- list(bp_mouse, go_mouse, kegg_mouse, reactome_mouse, wiki_pathways_mouse)
        return(de_enrichment_calc)
      }else if(input$selection_organism == "Zebrafish"){
        
        load("enrichment_sets/zebrafish.RData")
        all_genes <- rownames(results_df()[[1]])
        
        selected_genes <- rownames(results_df()[[1]][input$pca_pseudo_bulk_results_table_rows_all, ])
        
        bp_zebrafish <- prepare_for_enrichment(bp_zebrafish$terms, bp_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        go_zebrafish <- prepare_for_enrichment(go_zebrafish$terms, go_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        kegg_zebrafish <- prepare_for_enrichment(kegg_zebrafish$terms, kegg_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        reactome_zebrafish <- prepare_for_enrichment(reactome_zebrafish$terms, reactome_zebrafish$mapping, all_genes, feature_name = "gene_symbol")
        wiki_pathways_zebrafish <- prepare_for_enrichment(wiki_pathways_zebrafish$terms, wiki_pathways_zebrafish$mapping, all_genes, feature_name = "text_label")
        
        bp_zebrafish <- functional_enrichment(all_genes, selected_genes, bp_zebrafish)
        go_zebrafish <- functional_enrichment(all_genes, selected_genes, go_zebrafish)
        kegg_zebrafish <- functional_enrichment(all_genes, selected_genes, kegg_zebrafish)
        reactome_zebrafish <- functional_enrichment(all_genes, selected_genes, reactome_zebrafish)
        wiki_pathways_zebrafish <- functional_enrichment(all_genes, selected_genes, wiki_pathways_zebrafish)
        de_enrichment_calc <- list(bp_zebrafish, go_zebrafish, kegg_zebrafish, reactome_zebrafish, wiki_pathways_zebrafish)
        return(de_enrichment_calc)
        
      }
      # else if(input$selection_organism == "Integrated"){
      #   load("enrichment_sets/human.RData")
      #   all_genes <- rownames(results_df()[[1]])
      #   
      #   selected_genes <- rownames(results_df()[[1]][input$pca_pseudo_bulk_results_table_rows_all, ])
      #   
      #   bp_human <- prepare_for_enrichment(bp_human$terms, bp_human$mapping, all_genes, feature_name = "gene_symbol")
      #   go_human <- prepare_for_enrichment(go_human$terms, go_human$mapping, all_genes, feature_name = "gene_symbol")
      #   kegg_human <- prepare_for_enrichment(kegg_human$terms, kegg_human$mapping, all_genes, feature_name = "gene_symbol")
      #   reactome_human <- prepare_for_enrichment(reactome_human$terms, reactome_human$mapping, all_genes, feature_name = "gene_symbol")
      #   wiki_pathways_human <- prepare_for_enrichment(wiki_pathways_human$terms, wiki_pathways_human$mapping, all_genes, feature_name = "gene_symbol")
      #   
      #   bp_human <- functional_enrichment(all_genes, selected_genes, bp_human)
      #   go_human <- functional_enrichment(all_genes, selected_genes, go_human)
      #   kegg_human <- functional_enrichment(all_genes, selected_genes, kegg_human)
      #   reactome_human <- functional_enrichment(all_genes, selected_genes, reactome_human)
      #   wiki_pathways_human <- functional_enrichment(all_genes, selected_genes, wiki_pathways_human)
      #   de_enrichment_calc <- list(bp_human, go_human, kegg_human, reactome_human, wiki_pathways_human)
      #   return(de_enrichment_calc)}
    })

    output$pseudo_enrichment_table <- renderDT({
      req(pseudo_enrichment_calc())
      if(input$pseudo_enrichment_type == "GO"){
        datatable(
          pseudo_enrichment_calc()[[2]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip' 
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(pseudo_enrichment_calc()[[2]], is.numeric)), digits = 3)
      }else if(input$pseudo_enrichment_type == "BP"){
        datatable(
          pseudo_enrichment_calc()[[1]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip' 
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(pseudo_enrichment_calc()[[1]], is.numeric)), digits = 3)
      }else if(input$pseudo_enrichment_type == "KEGG"){
        datatable(
          pseudo_enrichment_calc()[[3]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip' 
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(pseudo_enrichment_calc()[[3]], is.numeric)), digits = 3)
      }else if(input$pseudo_enrichment_type == "Reactome"){
        datatable(
          pseudo_enrichment_calc()[[4]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip' 
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(pseudo_enrichment_calc()[[4]], is.numeric)), digits = 3)
      }else if(input$pseudo_enrichment_type == "WikiPathways"){
        datatable(
          pseudo_enrichment_calc()[[5]], 
          options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'l<"toolbar">rtip' 
          ),
          filter = 'top',
          rownames = FALSE,
          selection = "single",
          extensions = 'Scroller'
        ) %>%
          formatRound(which(sapply(pseudo_enrichment_calc()[[5]], is.numeric)), digits = 3)
      }
    })
    
    output$pseudo_ora_selection_ui <- renderUI({
      req(results_df())
      fluidRow(column(width = 6,
                      selectInput("pseudo_ora_selection", label = NULL, choices = c("CollecTRI", "PROGENy", "MSigDB Hallmark"))),
               column(width = 6,
                      actionButton("pseudo_run_ora", class = "btn-primary", "Calculate", width = '100%')))

    })
    
    collectri <- eventReactive(input$pseudo_run_ora,{
      req(input$pseudo_ora_selection == "CollecTRI")
      collectri <- read_rds("enrichment_sets/collectri.rds")
      tf_acts <- dc$run_ulm(mat = results_df()[[2]], net = collectri)
      tf_pvals <- tf_acts[[2]]
      tf_acts <- tf_acts[[1]]
      rownames(tf_pvals) <- "Selected_Clusters"
      rownames(tf_acts) <- "Selected_Clusters"
      results <- list(tf_acts, tf_pvals, collectri)
      return(results)
    })
    
    output$collectri_outputs <- renderUI({
      req(collectri())
      fluidRow(column(width = 6,
                      selectInput("collectri_plot_selection", label = NULL, choices = c("Top 25", "Volcano", "Network"))),
               column(width = 6,
                      conditionalPanel(condition = "input.collectri_plot_selection == 'Volcano'",
                                       selectizeInput("collectri_volcano_name", label = NULL, choices = colnames(collectri()[[1]]))),
                      conditionalPanel(condition = "input.collectri_plot_selection == 'Network'",
                                       selectizeInput("collectri_network_names", label = NULL, choices = colnames(collectri()[[1]]), multiple =T)
                                       )))
    })
    
    output$progeny_outputs <- renderUI({
      req(progeny())
      fluidRow(column(width = 6,
                      selectInput("progeny_plot_selection", label =NULL, choices = c("Scores", "Targets"))),
               column(width = 6,
                      conditionalPanel(condition = "input.progeny_plot_selection == 'Targets'",
                                       selectizeInput("progeny_scatter_selection", label = NULL, choices = c("Androgen", "EGFR", "Estrogen", "Hypoxia", "JAK-STAT", "MAPK",
                                                                                                           "NFkB", "PI3K", "TGFb", "TNFa", "Trail", "VEGF", "WNT", "p53"))
                                       )))
    })
    
    output$msigdb_outputs <- renderUI({
      req(msigdb())
      fluidRow(column(width = 6,
                      selectInput("msigdb_plot_selection", label =NULL, choices = c("Top Signatures", "Running Score"))),
               column(width = 6,
                      conditionalPanel(condition = "input.msigdb_plot_selection == 'Running Score'",
                                       selectizeInput("msigdb_running_selection", label = NULL, choices = sort(unique(msigdb()[[2]]$geneset)))
                      )))
    })
    
    
    output$collectri_top <- renderImage({
      req(collectri())
      dc$plot_barplot(
        acts=collectri()[[1]],
        contrast='Selected_Clusters',
        top= as.integer(25),
        vertical=F,
        # figsize=c(6, 10),
        save = "collectri.png"
      )
      list(src = "collectri.png")
    }, deleteFile = TRUE)
    
    output$collectri_volcano <- renderImage({
      req(collectri())
      logFCs = as.data.frame(t(results_df()[[1]]['log2FoldChange']))
      pvals = as.data.frame(t(results_df()[[1]]['padj']))
      rownames(logFCs) <- "Selected_Clusters"
      rownames(pvals) <- "Selected_Clusters"
      
      dc$plot_volcano(
        logFCs=logFCs,
        pvals=pvals,
        contrast='Selected_Clusters',
        name=input$collectri_volcano_name,
        net=collectri()[[3]],
        top=as.integer(10),
        sign_thr=0.05,
        lFCs_thr=0.5,
        save = "collectri_volcano.png"
      )
      list(src = "collectri_volcano.png")
    }, deleteFile = TRUE)
    
    
    output$collectri_network <- renderImage({
      req(collectri(),!is.null(input$collectri_network_names))
      dc$plot_network(
        net=collectri()[[3]],
        obs=results_df()[[2]],
        act=collectri()[[1]],
        n_sources= input$collectri_network_names,
        n_targets= as.integer(15),
        node_size=as.integer(100),
        figsize = c(5,5),
        c_pos_w='darkgreen',
        c_neg_w='darkred',
        vcenter=T,
        save = "collectri_network.png"
      )

      list(src = "collectri_network.png")
    }, deleteFile = TRUE)
    
    progeny <- eventReactive(input$pseudo_run_ora,{
      req(input$pseudo_ora_selection == "PROGENy")
      progeny <- read_rds("enrichment_sets/progeny.rds")
      mat = results_df()[[2]]
      
      pathway_acts = dc$run_mlm(mat=mat, net=progeny)
      pathway_acts <- pathway_acts[[1]]
      rownames(pathway_acts) <- "Selected_Clusters"
      results <- list(pathway_acts, progeny)
      return(results)
    })
    
    output$progeny_top <- renderImage({
      req(progeny())
      dc$plot_barplot(
        acts=progeny()[[1]],
        contrast='Selected_Clusters',
        top=as.integer(25),
        vertical=F,
        save = "progeny.png"
      )
      
      list(src = "progeny.png")
    }, deleteFile = TRUE)
    
    output$progeny_targets <- renderImage({
      req(progeny())
      dc$plot_targets(
        data=results_df()[[1]],
        stat='stat',
        source_name= input$progeny_scatter_selection,
        net=progeny()[[2]],
        top=as.integer(15),
        save = "progeny_targets.png"
      )
      
      list(src = "progeny_targets.png")
    }, deleteFile = TRUE)
    
    msigdb <- eventReactive(input$pseudo_run_ora,{
      req(input$pseudo_ora_selection == "MSigDB Hallmark")
      
      top_genes = results_df()[[1]][results_df()[[1]]['padj'] < 0.05,]
      msigdb <- read_rds("enrichment_sets/msigdb.rds")
      enr_pvals = dc$get_ora_df(
        df=top_genes,
        net=msigdb,
        source='geneset',
        target='genesymbol'
      )
      
      enr_pvals <- enr_pvals[order(enr_pvals$`Combined score`,decreasing = T),]
      results <- list(enr_pvals, msigdb)
      saveRDS(enr_pvals,"asdasds.rds")
      return(results)
    })
    
    output$msigdb_top <- renderImage({
      req(msigdb())
      
      dc$plot_dotplot(
        head(msigdb()[[1]], 15),
        x='Combined score',
        y='Term',
        s='Odds ratio',
        c='FDR p-value',
        scale=0.1,
        save = "msigdb_top.png"
      )
      
      list(src = "msigdb_top.png")
    }, deleteFile = TRUE)
    
    
    output$msigdb_running <- renderImage({
      req(msigdb(), input$msigdb_running_selection)
      
      dc$plot_running_score(
        df=results_df()[[1]],
        stat='stat',
        net=msigdb()[[2]],
        source='geneset',
        target='genesymbol',
        set_name= input$msigdb_running_selection,
        save = "msigdb_running.png"
      )
      
      list(src = "msigdb_running.png")
    }, deleteFile = TRUE)
    
    
    
    
    
    
    
}

options(shiny.sanitize.errors = TRUE,shiny.port=6868,shiny.host='0.0.0.0')

shinyApp(ui = ui, server = server)
