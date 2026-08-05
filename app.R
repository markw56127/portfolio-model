# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Author: Stephen Smith

library(shiny)
#library(pdfetch)
library(quantmod)
library(DT)
library(jsonlite)
library(httr)
library(tidyverse)

.pdenv <- new.env(parent=emptyenv())

pdfetch_new <- function(identifiers,
                        fields=c("open","high","low","close","adjclose","volume"),
                        from=as.Date("2007-01-01"),
                        to=Sys.Date(),
                        interval="1d") {
  
  valid.fields <- c("open","high","low","close","adjclose","volume")
  interval <- match.arg(interval, c("1d","1wk","1mo"))
  
  if (!missing(from))
    from <- as.Date(from)
  if (!missing(to))
    to <- as.Date(to)
  
  if (missing(fields))
    fields <- valid.fields
  if (length(setdiff(fields,valid.fields)) > 0)
    stop(paste0("Invalid fields, must be one of ", valid.fields))
  
  results <- list()
  from <- as.numeric(as.POSIXct(from))
  to <- as.numeric(as.POSIXct(to))
  
  # The following .yahooSession function is directly copied from quantmod, thank you to Joshua Ulrich.
  
  .yahooSession <- function(is.retry = FALSE) {
    cache.name <- "_yahoo_curl_session_"
    ses <- get0(cache.name, .pdenv) # get cached session
    
    if (is.null(ses) || is.retry) {
      ses <- list()
      ses$h <- curl::new_handle()
      # yahoo finance doesn't seem to set cookies without these headers
      # and the cookies are needed to get the crumb
      curl::handle_setheaders(ses$h, 
                              accept = "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7",
                              "User-Agent" = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36 Edg/115.0.1901.183")
      URL <- "https://finance.yahoo.com/"
      r <- curl::curl_fetch_memory(URL, handle = ses$h)
      # yahoo redirects to a consent form w/ a single cookie for GDPR:
      # detecting the redirect seems very brittle as its sensitive to the trailing "/"
      ses$can.crumb <- ((r$status_code == 200) && (URL == r$url) && (NROW(curl::handle_cookies(ses$h)) > 1))
      assign(cache.name, ses, .pdenv) # cache session
    }
    
    if (ses$can.crumb) {
      # get a crumb so that downstream callers don't have to handle invalid sessions.
      # this is a network hop, but very lightweight payload
      n <- if (unclass(Sys.time()) %% 1L >= 0.5) 1L else 2L
      query.srv <- paste0("https://query", n, ".finance.yahoo.com/v1/test/getcrumb")
      r <- curl::curl_fetch_memory(query.srv, handle = ses$h)
      if ((r$status_code == 200) && (length(r$content) > 0)) {
        ses$crumb <- rawToChar(r$content)
      } else {
        # we were unable to get a crumb
        if (is.retry) {
          # we already did a retry and still couldn't get a crumb with a new session
          stop("unable to get yahoo crumb")
        } else {
          # we tried to re-use a session but couldn't get a crumb
          # try to get a crumb using a new session
          ses <- .yahooSession(TRUE)
        }
      }
    }
    
    return(ses)
  }
  
  session <- .yahooSession()
  
  results <- data.frame()
  ### Added
  for (i in 1:length(identifiers)) {
    url <- paste0("https://query2.finance.yahoo.com/v8/finance/chart/",identifiers[i],
                  "?period1=",from,"&period2=",to,"&interval=",interval,"&events=history&crumb=",
                  session$crumb)
    response <- GET(url)
    list_response_text <- content(response,"text") 
    flat_text <- fromJSON(list_response_text,flatten = T)
    Result <- flat_text$chart$result
    
    DF <- data.frame(Date = as.character(as.Date.POSIXct(Result$timestamp[[1]])),
                     SYMBOL = Result$meta.symbol,
                     ADJ_CLOSE = round(Result$indicators.adjclose[[1]]$adjclose[[1]],2))
    
    results <- rbind(results, DF)
  }
  
  results_new <- results %>% 
    pivot_wider(names_from = SYMBOL,
                          values_from = ADJ_CLOSE)   
  results_new
}

ui <- fluidPage(
   
   # Application title
   titlePanel("Stock Portfolio"),
   p("Created by Stephen Smith and DuoDuo Ying"),
   p("PhD Students, UCLA Statistics"),
   p("Data for this page is gathered from"),
   a(href="https://finance.yahoo.com/","Yahoo Finance"),
   sidebarLayout(

      sidebarPanel(
        dateRangeInput('dateRange',
                       label = 'Date range input:',
                       start = Sys.Date() - 60, end = Sys.Date(),
                       format = "yyyy-mm-dd"
        ),
        textInput('stockChoices','Stock Choices (separate with commas)','^GSPC'),
        selectInput('frequency',
                    'Frequency',
                    choices = c('Daily','Weekly','Monthly'),
                    selected = 'Monthly'),
        submitButton('Get Adjusted Close Prices'),
        downloadButton("download","Download Data Table"),
        
        textOutput('unaccepted'),
        textOutput('badtickers')
      ),
      
      
      
      mainPanel(
         dataTableOutput('table')
      )
   )
)


server <- function(input, output) {

  unaccepted <- reactive({
    
    stocks <- input$stockChoices
    stocks <- gsub("\\s", "", stocks)
    choices <- unlist(strsplit(stocks,","))
    
    int <- input$frequency
    if (int == "Daily"){
      int <- "1d"
    }else if (int == "Weekly"){
      int <- "1wk"
    } else int <- "1mo"
    
    n <- nrow(
      na.omit(
        pdfetch_new('^GSPC',
                      from = input$dateRange[1],
                      to=input$dateRange[2],
                      interval = int,
                      fields = 'adjclose')))
    unaccepted <- c()
    badtickers <- c()
    accepted <- c()
    for (i in choices) {
      if (i == "") next
      df <- pdfetch_new(i,
                          from = input$dateRange[1],
                          to=input$dateRange[2],
                          interval = int,
                          fields = 'adjclose')

      if (is.null(df)){
        badtickers <- c(badtickers,i)
      } else if (nrow(na.omit(df))!=n) { # Inappropriate submission
        unaccepted <- c(unaccepted,i)
      } else accepted <- c(accepted,i)
    }
    return(list("unaccepted"=unaccepted,"badtickers"=badtickers,"accepted"=accepted))
  })
  
  stockInput <- reactive({
    # Save appropriate tickers
    choices <- unaccepted()[["accepted"]]
    
    stockdf <- data.frame()
    int <- input$frequency
    if (int == "Daily"){
      int <- "1d"
    }else if (int == "Weekly"){
      int <- "1wk"
    } else int <- "1mo"
    
    stockdf <- data.frame()
    for (i in choices) {
      df <- pdfetch_new(i,
                          from = input$dateRange[1],
                          to=input$dateRange[2],
                          interval = int,
                          fields = 'adjclose')
      df <- na.omit(df)
      
      if (i==choices[1]) stockdf <- data.frame("Date" = df[,1])
      stockdf <- cbind(stockdf, df[,2])
    }
    
    return(stockdf)
  })
  output$table <- DT::renderDataTable(stockInput())
  output$download <- downloadHandler(
    filename = "stockData.csv",
    content = function(file){
      write.csv(stockInput(),file)
    }
  )
  
  output$unaccepted <- renderText({
    tickers <- unaccepted()[[1]]
    
    if (length(tickers)>1) {
      val <- "There was insufficient data for "
      val2 <- paste(tickers,collapse = ", ")
      val3 <- ". These stocks were dropped from the final table."
      return(paste0(val,val2,val3))
    } else if (length(tickers)==1) {
      return(paste0("\nThere was insufficient data for\n",
                    tickers,
                    "\nThis stock was dropped from the final table."))
    } else return("")
  })
  
  output$badtickers <- renderText({
    tickers <- unaccepted()[[2]]
    
    if (length(tickers)>1) {
      val <- "The following tickers were not found: "
      val2 <- paste(tickers,collapse = ", ")
      val3 <- "."
      return(paste0(val,val2,val3))
    } else if (length(tickers)==1) {
      return(paste0("The following ticker was not found: ",tickers,"."))
    } else return("")
  })
  

}

# Run the application 
shinyApp(ui = ui, server = server)

