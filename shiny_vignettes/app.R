library(shiny)
library(stringr)

# what is shown on screen, and the user can interact with.
ui <- fluidPage(
  # Head with selectors and link to dataset
  
  fluidRow(
    tabsetPanel(
      tabPanel("Info", style= "margin-left: 10px;
               margin-right: 10px;",
               htmlOutput("info")
      ),
      tabPanel("vignettes",
               fluidRow(
                 column(3),
                 column(3, align = "center", 
                        selectInput("vignette", label = "Select Vignette to view", choices = list(""), selected = NULL, width = "400px"),
                 ),
                 column(3, style = "margin-top: 25px;",
                        column(6, align="center",
                               downloadButton("downloadHtml", "Download HTML")),
                        column(6, align = "center",
                               downloadButton("downloadRmd", "Download original file"))),
                 column(3)),
               fluidRow(style= "margin-left: 10px;
                        margin-right: 10px;",
                 htmlOutput("inc"))
      )
  ))
)

# All the plotting and loading processes
server <- function(input, output, session) {
  
  # get all rds files in vignettes dir
  # here, need files to be in vignettes/name_of_directory/ and in format '.html'
  # exemple : ./vignettes/cellranger/cellranger_tutorial_neuron1k.nb.html
  # No vignette should have the same name, even if in different directory.
  # step ran only at start
  allHtmlFiles <- reactive({
    list.files('./vignettes/', pattern = '.html', include.dirs = T, recursive = T)
  })
  
  allRmdFiles <- reactive({
    list.files('./vignettes/', pattern = '.[Rr]md|.ipynb', include.dirs = T, recursive = T)
  })
  
  # parse name of subdir + file to get a named list like list$subdir = c(file1, file2)
  # step ran only at start
  inputNamed <- reactive({
    x <- allHtmlFiles()
    # data.table col1 dir col2 file
    df <- str_split_fixed(x, pattern = "/", n = str_count(x, "/") + 1)
    # strip filename from ".rds"
    df[,2] <- gsub(".nb.html$|.html$", '', df[,2])
    df[,2] <- gsub("_", ' ', df[,2])
    df[,1] <- gsub("_", ' ', df[,1])
    # make a named lists, names being dir and items files 
    x <- split(df[,2], df[,1])
    # Even with just one item, make all items as list
    for (v in names(x)){
      x[[v]] = as.list(x[[v]])
    }
    x
  })
  
  observe({
    updateSelectInput(session, "vignette", label = "Select Vignette to view", choices = c("", inputNamed()), selected = "")
  })
  
  getHtmlFile <- reactive({
    validate(need(input$vignette, message = "Please select a valid vignette"))
    fileName = paste('./vignettes/', grep(gsub(" ", "_", input$vignette[1]), allHtmlFiles(), value = T), sep = "")
    return(fileName)
  })
  
  getRmdFile <- reactive({
    validate(need(input$vignette, message = "Please select a valid vignette"))
    fileName = paste('./vignettes/', grep(gsub(" ", "_", input$vignette[1]), allRmdFiles(), value = T), sep = "")
    return(fileName)
  })
  
  getPage<-function() {
    return(includeHTML(getHtmlFile()))
  }
  
  output$inc<-renderUI({getPage()})
  
  getPresentation<-function() {
    return(includeHTML("./presentation.html"))
  }
  
  output$info<-renderUI({getPresentation()})
  
  output$downloadRmd <- downloadHandler(
    filename <- function(){
      paste0(gsub(" ", "_", input$vignette[1]), "_out", str_extract(getRmdFile(), ".*"))
    },
    content <- function(file) {
      file.copy(getRmdFile(), file)
    }
  )
  
  output$downloadHtml <- downloadHandler(
    filename <- function(){
      paste0(gsub(" ", "_", input$vignette[1]), "_out", str_extract(getHtmlFile(), ".*"))
    },
    content <- function(file) {
      file.copy(getHtmlFile(), file)
    }
  )

}

shinyApp(ui = ui, server = server)