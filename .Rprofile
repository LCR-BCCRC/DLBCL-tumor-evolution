source("renv/activate.R")

# Added per https://github.com/REditorSupport/vscode-R/wiki/R-Session-watcher#advanced-usage-for-self-managed-r-sessions
# Trying to get Radian to use the VScode viewer instead of X11
Sys.setenv(TERM_PROGRAM = "vscode")

if (interactive() && Sys.getenv("RSTUDIO") == "") {
  source(file.path(Sys.getenv(if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"), ".vscode-R", "init.R"))
}

# Watch global environemnt symbols to provide hover on and completion after session symbol.
# Only specify in .Rprofile since it only takes effect on session startup.
options(vsc.globalenv = TRUE)

# Which view column to show the plot file on graphics update
# Use FALSE to diable plot watcher so that the default R plot device is used.
# Only specify in .Rprofile since it only takes effect on session startup.
if ("httpgd" %in% .packages(all.available = TRUE)) {
  options(vsc.plot = FALSE)
  options(device = function(...) {
    httpgd::httpgd()
    .vsc.browser(httpgd::httpgdURL(), viewer = "Beside")
  })
} else {
    options(vsc.plot = "Beside")
}


# The arguments for the png device to replay user graphics to show in VSCode.
# Ignored if options(vsc.plot = FALSE).
options(vsc.dev.args = list(width = 800, height = 600))

# Which view column to show the WebView triggered by browser (e.g. shiny apps)?
# Use FALSE to open in external web browser.
options(vsc.browser = "Active")

# Which view column to show the WebView triggered by viewer (e.g. htmlwidgets)?
# Use FALSE to open in external web browser.
options(vsc.viewer = "Two")

# Which view column to show the WebView triggered by page_viewer (e.g. profvis)?
# Use FALSE to open in external web browser.
options(vsc.page_viewer = "Active")

# Which view column to show the WebView triggered by View()?
# Use FALSE for R's native View(), which should be specified in .Rprofile
#   since it only takes effect on session startup.
options(vsc.view = "Two")

# Which view column to show the WebView triggered by help panel
# (e.g. after sending `?mean` to terminal)?
# Use FALSE to disable help panel and revert to old behaviour.
options(vsc.helpPanel = "Two")

# How much of the object to show on hover and autocomplete detail?
# As controlled by max.level arg of str().
# Use 0 (or 1) is the default - literal value or object type and dimensions
# Use 2 to show list contents, data frame columns, and example values.
options(vsc.str.max.level = 0)
