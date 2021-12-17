source_all = function(dir_name){
	
	files = dir(dir_name)
	for(f in 1:length(files)){
		source(here(dir_name, files[f]))
	}
	
}