data_manipulation <- function(original_data, original_model, segments) {
  # Creating variables that contain the names of measurement items and constructs
  item_names <- seminr:::all_items(original_model$measurement_model)
  construct_names <- seminr:::construct_names(original_model$structural_model)
  
  # Extracting viable segment IDs
  viable_segments_ID <- as.character(segments$viable) 
  
  # Creating a list containing the node cases of the viable segments
  viable_segments <- sapply(viable_segments_ID,function(id) node_cases = segments$node_cases[[id]])
  
  # Dividing the original pooled data into separate datasets for each segment
  segment_dataset <- function(dataset, segments) {
    lapply(segments, function(segment) {
      dataset[segment, ]
    })
  }
  
  segmented_datasets <- segment_dataset(original_data, viable_segments)
  
  # Creating a vector with group names
  segment_labels <- rep(NA, nrow(original_data)) #create a new vector to store the segment_id
  
  # Add the value of segment_id to each observation (row)
  sapply(seq_along(viable_segments), function(i) {
    segment_labels[as.numeric(viable_segments[[i]])] <<- paste0("segment_", viable_segments_ID[i])
  })
  
  return(list(item_names = item_names, construct_names = construct_names, viable_segments_ID = viable_segments_ID, segment_labels = segment_labels))
}