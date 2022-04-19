# # Run once
# library(remotes)
# install.packages("extrafont")
# install.packages("cairo_pdf")
# remotes::install_version("Rttf2pt1", version = "1.3.8")
# library(extrafont)
# # Install **TTF** Latin Modern Roman fonts from www.fontsquirrel.com/fonts/latin-modern-roman
# # Import the newly installed LModern fonts, change the pattern according to the
# # filename of the lmodern ttf files in your fonts folder
# font_import(pattern = "lmroman*")

# Run each time
library(extrafont)
loadfonts(device = "win")

## Create an epiweek calendar
dates <- seq(as.Date("2020-01-01"),as.Date("2021-04-01"),by="1 day")
date_start <- rep(as.Date("2020-01-01"), times=length(dates))
epiweeks <- floor(difftime(dates,date_start,units = "weeks"))+1
epi_calendar <- tibble(date=dates,week=epiweeks)
epi_calendar <- epi_calendar %>% group_by(week) %>% mutate(first_day=min(date))

export_theme <- theme_tufte() +
  theme(
    axis.text.x = element_text(size=7,family="LM Roman 10"),
    axis.text.y=element_text(size=7,family="LM Roman 10"),
    axis.title.x=element_text(size=10,family="LM Roman 10",vjust=-1),
    axis.title.y=element_text(size=10,family="LM Roman 10"),
    
    ## Axis lines
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),
    
    ## Title
    plot.title = element_text(family="LM Roman 10",size=8,face="bold",hjust=0.5),
    plot.tag = element_text(family="LM Roman 10",size=10,face="bold"),
    
    ## Legends
    legend.title=element_text(size=8,family="LM Roman 10",face="italic"),
    legend.text=element_text(size=8,family="LM Roman 10"),
    legend.key.size= unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=8,family="LM Roman 10",face="bold"),
    #strip.background=element_rect(fill="#f0f0f0")
    strip.background=element_blank()
  )



AAAS_palette <- c("blue1"="#3B4992FF","red1"="#EE0000FF","green1"="#008B45FF",
                  "purple1"="#631879FF","teal1"="#008280FF","red2"="#BB0021FF",
                  "purple2"="#5F559BFF","purple3"="#A20056FF",
                  "grey1"="#808180FF","black"="#1B1919FF")