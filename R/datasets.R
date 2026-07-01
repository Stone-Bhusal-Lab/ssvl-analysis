add_dataset <- function(store, dataset) {
  
  current <- store()
  
  current[[dataset$name]] <- dataset
  
  store(current)
}

remove_dataset <- function(store, name) {
  
  current <- store()
  
  current[[name]] <- NULL
  
  store(current)
}

rename_dataset <- function(store, old, new) {
  
  current <- store()
  
  current[[new]] <- current[[old]]
  
  current[[old]] <- NULL
  
  current[[new]]$name <- new
  
  store(current)
}