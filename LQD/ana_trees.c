void ana_trees(int jobnumber){
  gROOT->LoadMacro("analyze_trees_c.so"); 
  
  analyze_trees(jobnumber);

}
