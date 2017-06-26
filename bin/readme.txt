Download the whole bin folder if you want to run the GUI file VSM_GUI_Win64.exe, which is depend on the dynamic link libraries.

Buttons on the top menu for constraints selection:
1. Open File
2. Switch the Main View and Detial View.
3. Node Dragger
4. Constraints on one node
  a) Equal Angle Constraint on the selected node if the degree of this node is no smaller than 3.
  b) One or two largest Circle Constraints on the selected node if there is circle(s) on this node.
  Note: First select a node, then click the constraint button.
5. Constraints on a lasso area
  a) Equal Angle Constraints on nodes which are in the lasso area.
  b) Horizontal Constraints on the edges whose both endpoints are in the lasso area.
  c) Vertical Constraints on the edges whose both endpoints are in the lasso area.
  d) User Assigned Direction Constraints on the edges whose both endpoints are in the lasso area.
  e) Crossing Removing Constraints on the edges whose both endpoints are in the lasso area.
  f) Orthogonal Constraints on the nodes which are in the lasso area and the degree is no bigger than 4.
  g) Symmetry Constraints on the lasso area.
  Note: For a), b), c), e), and f), first, draw a lasso area; Then, click the constraint button.
        For d) and g), first, draw lasso area; Then, click Constraint button; Finally, you need to draw a axis as direction-axis or symmetry-axis.     
6. Community(Cluster) Non-overlap Constraint
  Note: Only for the data with community information.
7. Constraints on two nodes
  a) Path-finder. 
     Direction constraints on the path in the Main View with idea direction sum(direction(i,j)), where i and j is adjacent nodes on the path
     Vertical constraints and Node Non-overlap Constraints on the path in the Detail View.
  b) Common-node finder
  Note: First, click the constraint button; Then, select source node and target node
 8. Reset
 
Detail View:
  Adjacent nodes will show up when user double click on one node in detial view. Corresponding nodes will be highlighted in the Main View.
  
Weight Constrolling View:
  Listing the constraints in the graph.
  User scale the weight of each constraint by the sliding block. 
  Apply the weight with the Confirm Button and delete the constraint with the Cancel Button.
  
