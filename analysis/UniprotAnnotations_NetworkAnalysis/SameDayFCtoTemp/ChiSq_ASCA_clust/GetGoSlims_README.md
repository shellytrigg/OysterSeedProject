### Getting GO Slim terms for a list of GO IDs

1. Make a list of your genes/proteins with GO IDs
	- merge your Uniprot BLAST output with a Uniprot master table (contains GO IDs)

	example BLAST output


	example Uniprot Master table

	example of merged tables

	- Subset only gene/protein IDs and GO IDs
	- Move all GO IDs into their own column
	- Reshape data so that all GO IDs are listed in one column with their corresponding gene/protein ID in the other column

2. Make a list of unique GO IDs 
3. Use GSEAbase to get GO slims for your GO IDs
	- convert your list of unique GO IDs to a GOCollection
	- load and format GO slim data for GSEAbase
	- get GO slims
