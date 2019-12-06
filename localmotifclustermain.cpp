// localmotifclustermain.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "localmotifcluster.h"
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <string>


//Returns most central node ID according to degree centrality in the transformed graph
int getMostCentralNode(ProcessedGraph* graph_p, std::map<int,bool>* coveredNodes){

  //Store transformed graph and weights of the graph
  PUNGraph transformedGraph = graph_p->getTransformedGraph();   
  WeightVH weights = graph_p->getWeights();
  
  //Keep current node ID and maxDegree which will be updated in one pass through all nodes
  int currNodeID;
  int maxDegree;
  maxDegree = -1;
  
  
  //Foreach node in the transformed graph, check if they do not exist in the coveredNodes map
  //If they dont, update maxDegree, currNodeID if necessary
  for (TUNGraph::TNodeI NI = transformedGraph->BegNI(); NI < transformedGraph->EndNI(); NI ++ ) {
    int nodeID = NI.GetId();
    int degree = int(weights[nodeID](nodeID));
    //If node not yet covered
    if (coveredNodes->find(nodeID) == coveredNodes->end()){
      //printf("NOT FOUND :(\n");
      if (degree > maxDegree){
	    currNodeID = nodeID;
	    maxDegree = degree;
	    //printf("Node %d with degree %d\n", nodeID, maxDegree);
	  }
    } 
  }
  return currNodeID;
}

//Reads the ground-truth community specified file
std::vector< std::map<int,bool> > getGroundTruthCommunities(std::string fname){
  std::vector< std::map<int,bool> > gtPartitioning;
  
  std::ifstream file(fname.c_str());
  
  if (file.is_open()){
    std::string line;
    while (getline(file, line)) {
	    if (line.empty())
      {
		    continue;
        //cout << "Empty line." << endl;
      }
	    std::string lineText = line.c_str();
      // using printf() in all tests for consistency
      //printf("%s\n", lineText);
      //std::cout << lineText << std::endl;
      std::string buf;                 // Have a buffer string
      std::stringstream ss(lineText);       // Insert the string into a stream
     
      std::map<int,bool> community; // Create vector to hold our words
      bool flag = false;
      while (ss >> buf){
        //tokens.push_back(buf.c_str());
        //if (!buf == " "){
	    //  std::cout << buf + "\n" << std::endl;
	    //}
        
		
        if(buf.find_first_not_of('\n') != std::string::npos)
        {
          // There's a non-space.
          //std::cout << buf << std::endl;
          std::string node(buf);
          int nodeID = atoi(buf.c_str());
          //std::cout << nodeID << std::endl;
          //Add nodeID to current GT community
          community.insert( std::pair<int,bool>(nodeID, true) );
          //Mark that the community is non-empty
          flag = true;
        }


      }
      if (flag){
        //Add GT community to global community partitioning
        gtPartitioning.push_back(community);
      }
      
    }
      
    file.close();
  }

  return gtPartitioning;
} 

//Compute Jaccard Similarity between two clusters c1,c2
double jaccSim(std::map<int,bool>* c1, std::map<int,bool>* c2){
  double intersectCounter = 0.0;
  double sizeC1 = double(c1->size()); double sizeC2 = double(c2->size());

  std::map<int,bool>* smallestC;
  std::map<int,bool>* largestC;

  //Loop over smallest of c1 and c2
  if (sizeC1 > sizeC2){
    smallestC = c2;
    largestC = c1;
  } else{
    smallestC = c1;
    largestC = c2;
  }
  //std::cout << "HERE NOW" << std::endl;
  for(std::map<int,bool>::iterator iter = smallestC->begin(); iter != smallestC->end(); ++iter)
  {
    
    int node =  iter->first;
    //std::cout << node << std::endl;
    if (largestC->find(node) != largestC->end()){
      //Found overlap
      intersectCounter += 1.0;
    }

    //ignore value
    //Value v = iter->second;
  }
  double unionCounter = sizeC1 + sizeC2 - intersectCounter;
  return intersectCounter/unionCounter;
}

//Computes and writes to memory precision, recall and f1 score
void evaluate(std::map<int,bool>* truth, std::map<int,bool>* pred, 
              double &precision, double &recall, double &f1){
  double TP = 0.0; double FP = 0.0; double FN = 0.0;
  double truthSize = double(truth->size()); double predSize = double(pred->size());
  //Loop over predicted community, and keep TP (present in truth) and FP (missing in truth) counter
  for(std::map<int,bool>::iterator iter = pred->begin(); iter != pred->end(); ++iter){
    int predNode = iter->first;
    if (truth->find(predNode) != truth->end()){
      //If node present in ground truth community, increment TP
      TP++;
    }else{
      //Node missing in ground truth community --> FP
      FP++;
    }
  }
  //Use the identity Len(truth) - TP = FN
  FN = truthSize - TP;

  precision = TP / (TP + FP);
  recall = TP / (TP + FN);
  if (precision == 0 && recall == 0){
    f1 = 0;
  }else{
    f1 = 2 * (precision * recall)/(precision + recall);
  }
  
}

//Write precision, recall, and f1 measure to file
void writeResults(double& precision, double& recall, double& f1, std::string resultFile
                  ,const char* runtime) {
  std::ofstream myfile;
  const char* writeTo = resultFile.c_str();
  myfile.open(writeTo);
  myfile << "Precision " << precision << "\n";
  myfile << "Recall " << recall << "\n";
  myfile << "F1 " << f1 << "\n";
  myfile << "Runtime " << runtime << "\n";
  myfile.close();
}

void writeRuntime(std::string runtimeFile, const char* runtime) {
	std::ofstream myfile;
	const char* writeTo = runtimeFile.c_str();
	myfile.open(writeTo);
	myfile << "Runtime " << runtime << "\n";
	myfile.close();
}

void writeClusters(std::vector< std::map<int,bool> >& clusters, std::string fout){
  std::ofstream myfile;
  const char* writeTo = fout.c_str();
  myfile.open(writeTo);
  for (std::vector< std::map<int,bool> >::iterator iter = clusters.begin(); iter != clusters.end(); iter++){
    std::map<int,bool> currClust = *iter;
    std::map<int, bool>::iterator nodeIter = currClust.begin();

    //Create line string 
    std::string lineToWrite = "";

    while (nodeIter != currClust.end()){
      int node = nodeIter->first;

      std::ostringstream s;
      s << "" << node;
      std::string nodeString(s.str());

      lineToWrite += nodeString;
      lineToWrite += " ";


      nodeIter++;
    }

    myfile << lineToWrite << "\n";

  } 
  myfile.close();

}

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("Local motif clustering. build: %s, %s. Time: %s",
       __TIME__, __DATE__, TExeTm::GetCurTm()));  
  TExeTm ExeTm;  
  Try

  const bool IsDirected = 
    Env.GetIfArgPrefixBool("-d:", false, "Directed graph?");

  const bool DanglingMerging = 
    Env.GetIfArgPrefixBool("-q:", false, "Dangling node merging?");

  //Ground-truth communities file
  std::string clustout_filename =
      Env.GetIfArgPrefixStr("-w:", "none.txt", "clustout").CStr();

  std::string runtime_filename =
	  Env.GetIfArgPrefixStr("-r:", "none.txt", "runtime").CStr();
  //const std::string fname = gt_filename.CStr();
  

  //printf("%s HERE IS YOUR ARGUMENT \n", gt_filename.CStr()); 
  //std::vector< std::map<int,bool> > gtCommunities = getGroundTruthCommunities(fname);



  ProcessedGraph graph_p;
  if (IsDirected) {
    const TStr graph_filename =
      Env.GetIfArgPrefixStr("-i:", "C-elegans-frontal.txt", "Input graph file");
    const TStr motif =
      Env.GetIfArgPrefixStr("-m:", "triad", "Motif type");
    MotifType mt = ParseMotifType(motif, IsDirected);
    PNGraph graph;
    if (graph_filename.GetFExt().GetLc() == ".ngraph") {
      TFIn FIn(graph_filename);
      graph = TNGraph::Load(FIn);
    } else if (graph_filename.GetFExt().GetLc() == ".ungraph") {
      TExcept::Throw("Warning: input graph is an undirected graph!!");
    } else {
      graph = TSnap::LoadEdgeList<PNGraph>(graph_filename, 0, 1);
    }
    TSnap::DelSelfEdges(graph);
    graph_p = ProcessedGraph(graph, mt);

  } else {

    const TStr graph_filename =
      Env.GetIfArgPrefixStr("-i:", "C-elegans-frontal.txt", "Input graph file");
    const TStr motif =
      Env.GetIfArgPrefixStr("-m:", "clique3", "Motif type");
    MotifType mt = ParseMotifType(motif, IsDirected);
    PUNGraph graph;
    if (graph_filename.GetFExt().GetLc() == ".ungraph") {
      TFIn FIn(graph_filename);
      graph = TUNGraph::Load(FIn);
    } else if (graph_filename.GetFExt().GetLc() == ".ngraph") {
      TExcept::Throw("Warning: input graph is a directed graph!!");
    } else {
      graph = TSnap::LoadEdgeList<PUNGraph>(graph_filename, 0, 1);
    }
    TSnap::DelSelfEdges(graph);
    graph_p = ProcessedGraph(graph, mt);
  }

  const TInt seed =
    Env.GetIfArgPrefixInt("-s:", 1, "Seed");
  const TFlt alpha =
    Env.GetIfArgPrefixFlt("-a:", 0.98, "alpha");
  const TFlt eps =
    Env.GetIfArgPrefixFlt("-e:", 0.0001, "eps");


  std::cout << "Graph loaded into memory" << std::endl;
  //Store number of nodes in graph
  int numOfNodes = graph_p.getOriginalGraph()->GetNodes();
  
  //printf("%d nodes\n",numOfNodes);
  
  //Map of covered nodes (key = nodeID, value = True (if covered))
  std::map<int,bool>  covered;
  //Stores clusters
  std::vector< std::map<int,bool> > partitioning;
  
  //use partitioning with     arr.push_back(cluster); for (auto x : partitioning) {...}
  
  //Do local clustering until entire graph is partitioned
  while (covered.size() != numOfNodes){
	//Get the most central node according to degree centrality (passing only pointers)
    int centralNode = getMostCentralNode(&graph_p, &covered);
    //printf("%d is the most central node\n", centralNode);  
	  
    MAPPR mappr;
    mappr.computeAPPR(graph_p, centralNode, alpha, eps / graph_p.getTotalVolume() * graph_p.getTransformedGraph()->GetNodes());
    //Sweep param: -1 first local min. 0: global min. 
    mappr.sweepAPPR(0);
    //mappr.printProfile();
    //printf("Size of Cluster: %d.\n", mappr.getCluster().Len());

    //The cluster found by MAPPR
    TIntV foundCluster = mappr.getCluster();
    std::map<int,bool> dictCluster;

    
    //Mark them as covered
    for (TInt* it = foundCluster.BegI(); it < foundCluster.EndI(); it++){
      //printf("%d node in cluster\n", int(*it));
      int nodeInt = int(*it);
      covered.insert(std::pair<int,bool>(nodeInt,true));
      dictCluster.insert(std::pair<int,bool>(nodeInt, true));
    }
    //Add cluster to partitioning
    partitioning.push_back(dictCluster);



  }
  
  std::cout << "Done partitioning" << std::endl;

  const char* runtime = ExeTm.GetTmStr();

  if (DanglingMerging){

	  /*
	  Should have used the for loop in line 274, but somehow it doesn't work.
	  Loop over every cluster, assign those with size 1 to other cluster.
	  The criteria is which cluster most of its neibors belong to.
	   */
	  int psize = partitioning.size();
	  //for (std::vector< std::map<int,bool> >::iterator vit = partitioning.begin(); vit != partitioning.end(); vit++){
	  for (int i = psize - 1; i >= 0; i--) {
		  std::map<int, bool> vre = partitioning[i];
		  //if (vit->size() == 1){
		  if (vre.size() == 1) {
			  std::map<std::map<int, bool>, int> candidate;
			  //int iKey = vit->begin()->first;
			  int iKey = vre.begin()->first;
			  TUNGraph::TNodeI NI = graph_p.getOriginalGraph()->GetNI(iKey);

			  /*
			  Check which cluster(s) its neibors belong to, and add them into candidate.
			  */
			  for (int j = 0; j < NI.GetOutDeg(); j++) {
				  int Nbr = NI.GetNbrNId(j);
				  for (std::vector< std::map<int, bool> >::iterator vit2 = partitioning.begin(); vit2 != partitioning.end(); vit2++) {
					  std::map<int, bool>  big_cluster = *vit2;
					  std::map<int, bool>::iterator it = big_cluster.find(Nbr);

					  /* the node is in this big cluster */
					  if (it != big_cluster.end()) {
						  if (candidate.find(big_cluster) != candidate.end()) {
							  candidate.find(big_cluster)->second += 1;
						  }
						  else {
							  candidate.insert(std::pair<std::map<int, bool>, int>(big_cluster, 1));
						  }
					  }
				  }
			  }
			  /* Someone only send an email to oneself, thus leave them alone. */
			  if (candidate.size() == 0) {
				  //cout << "node:" << iKey << " degree: " << NI.GetDeg() << endl;
				  continue;
			  }

			  /* determine the final cluster from candidates */
			  int maxOcurrence = 0;
			  std::map<std::map<int, bool>, int>::iterator maxKey;

			  for (std::map<std::map<int, bool>, int>::iterator it = candidate.begin(); it != candidate.end(); it++) {
				  if (it->second > maxOcurrence)
					  maxKey = it;
				  maxOcurrence = it->second;
			  }

			  std::vector< std::map<int, bool> >::iterator desiredCls = find(partitioning.begin(), partitioning.end(), maxKey->first);
			  if (desiredCls != partitioning.end()) {
				  /* add the node into the final disired cluster */
				  (*desiredCls).insert(std::pair<int, bool>(iKey, true));

				  /* delete the cluster with that node from partitioning */
				  std::map<int, bool>discardMap;
				  discardMap.insert( std::pair<int,bool>(iKey, true) );
				  std::vector< std::map<int, bool> >::iterator discardCls = find(partitioning.begin(), partitioning.end(), discardMap);
				  if (discardCls != partitioning.end()) {
					  partitioning.erase(discardCls);
				  }
				  else {
					  std::cout << "something wrong 1" << std::endl;
				  }
			  }
			  else {
				  std::cout << "something wrong 2" << std::endl;
			  }
		  }
	  }




  /*for (int i = 0; i < partitioning.size(); i++) {
	  if (partitioning[i].size() == 1) {
		  cout << "size of 1: " << partitioning[i].begin()->first << endl;
	  }
  }*/
  }



  //TODO:merging candidate clusters?... nah maybe later


  /*
  double macroPrec = 0.0, macroRec = 0.0, macroF1 = 0.0;

  //TODO: Compare candidate clusters to GT clusters 
  for(int i = 0; i < partitioning.size(); i++){
	  std::map<int,bool> outCluster = partitioning[i];

    //Max jacc sim so far
    double maxSim = -1;
    //Best ground cluster (max overlap) so far
    std::map<int,bool>* groundCluster;

    //Loop over each GT community and find one with largest overlap
    for (int j = 0; j < gtCommunities.size(); j++){
      std::map<int,bool> currGroundCluster = gtCommunities[j];
      double sim = jaccSim(&outCluster, &currGroundCluster);
      if (sim > maxSim){
        maxSim = sim;
        groundCluster = &gtCommunities[j];
        //std::cout << "Updated cluster to " + j << std::endl;
      }
    }

    //Compute Precision, Recall, F1 for outCluster and groundCluster
    
    double precision, recall, f1;
    evaluate(groundCluster,&outCluster,precision, recall, f1);
    /*std::cout << "Precision" << std::endl;
    std::cout << precision << std::endl;
    std::cout << "Recall" << std::endl;
    std::cout << recall << std::endl;
    std::cout << "F1" << std::endl;
    std::cout << f1 << std::endl;
    macroPrec += precision;
    macroRec += recall;
    macroF1 += f1;
  }
  macroPrec /= double(partitioning.size());
  macroRec /= double(partitioning.size());
  macroF1 /= double(partitioning.size());
  
  std::cout << "Precision" << std::endl;
  std::cout << macroPrec << std::endl;
  std::cout << "Recall" << std::endl;
  std::cout << macroRec << std::endl;
  std::cout << "F1" << std::endl;
  std::cout << macroF1 << std::endl;
  

  
  writeResults(macroPrec, macroRec, macroF1, write_filename, runtime);
  */
  writeClusters(partitioning, clustout_filename);
  writeRuntime(runtime_filename, runtime);

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(),
   TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
