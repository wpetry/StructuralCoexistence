library(shiny)
library(shinyIncubator)
library(rgl)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    
    # Application title
    titlePanel("Understanding Structural Coexistence for Multispecies Communities"),
    
    # Inputs
    
    sidebarLayout(
     
      sidebarPanel(
        radioButtons("spp", "Number of Species",
                     choices = c('3','4'),
                     selected = '3'),
        conditionalPanel(
          condition = "input.spp == '3'",
          numericInput("r1","r1",1,min=0,max=1,step=0.1),
          numericInput("r2","r2",1,min=0,max=1,step=0.1),
          numericInput("r3","r3",0.5,min=0,max=1,step=0.1),
          matrixInput("alphamat",
                      "Competition coefficients matrix (α)",
                      data.frame(matrix(c(1,0.5,0.05,0.4,1,0.5,0.3,0.6,1),nrow=3))),
          checkboxInput("draw_table","Show table",value=T),
          radioButtons("cone_opt", "Plot cone",
                       choices = c("none", "static", "interactive"),
                       selected = "none")
        ),
        conditionalPanel(
          condition = "input.spp == '4'",
          numericInput("rr1","r1",1,min=0,max=1,step=0.1),
          numericInput("rr2","r2",1,min=0,max=1,step=0.1),
          numericInput("rr3","r3",0.5,min=0,max=1,step=0.1),
          numericInput("rr4","r4",1,min=0,max=1,step=0.1),
          matrixInput("alphamat4",
                      "Competition coefficients matrix (α)",
                      data.frame(matrix(c(1,0.25,0.05,0.1,0.2,1,0.2,0.3,0.3,0.05,1,0.3,0.3,0.2,0.4,1),nrow=4))),
          checkboxInput("draw_table4sp","Show table",value=T)
          
        )
   
       ),
    
      # Outputs
      
     mainPanel(
        conditionalPanel("input.spp == '3' & input.draw_table",
                         fluidRow(tableOutput("stats")),
                         align = "center"),
        conditionalPanel("input.spp == '3'",
                         plotOutput("proj"),
                         align = "center"),
        conditionalPanel("input.spp == '3' & input.cone_opt == 'static'",
                         plotOutput("cone"), align = "center"),
        conditionalPanel("input.spp == '3' & input.cone_opt == 'interactive'",
                         rglwidgetOutput("cone3d"), align = "center"),
        conditionalPanel("input.spp == '4' & input.draw_table4sp",
                         fluidRow(tableOutput("stats4sp")),
                         align = "center"),
        conditionalPanel("input.spp == '4'",
                         rglwidgetOutput("proj4sp"), align = "center")
      )
     
    )
  ))

