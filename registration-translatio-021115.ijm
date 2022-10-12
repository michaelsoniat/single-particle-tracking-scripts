//   Ajust contrast of all frames. 
// The output movie has all the pixel values changed permanently. 
// Output movie is saved on the same folder as original movie. 
// Keep all your original files. 

pathref = File.openDialog("Open the fix point text file");  // opens a text file for reference translations. be sure this fixed point blinks as little as possible.
filestring=File.openAsString(pathref); 
rows=split(filestring, "\n"); 
cref=newArray(rows.length); 
xref=newArray(rows.length); 
yref=newArray(rows.length); 
for(i=0; i<rows.length; i++){ 
columns=split(rows[i],"\t"); 
cref[i]=parseFloat(columns[0]); 
xref[i]=parseFloat(columns[1]); 
yref[i]=parseFloat(columns[2]); 

}   // end of text file opening

 path1 = File.openDialog("Select a File");  // this opens the file you want to translate according to the fixed point you selected before. 
dir = File.getParent(path1);
open(path1);
n=nSlices;
Array.getStatistics(cref, min, max, mean, stdDev);
print(max);
// Create a directory in temp
  myDir = dir+File.separator+"my-test-dir"+File.separator;
  File.makeDirectory(myDir);            // creates a temporary file where all frames will be processed independently. 
  if (!File.exists(myDir))
      exit("Unable to create directory");
  print("");
 

print("Original movie:");
print(path1);


run("Image Sequence... ", "format=TIFF name=C2-test-20 start=0 digits=4 save="+myDir);  // this converts the stack into frames inside myDir

close();
  
requires("1.39l");
 // nDirs =1;
  startingImage = 1;
  increment = 1;
   dir = myDir;    
destination=dir;



  index = 0;
//sa=0;
xorigin=xref[0];    // all frames in the movie are translated back at the origin. This is the origin. 
yorigin=yref[0];
inpast=0;



list = getFileList(dir);   // each file is a frame of the movie. 

      for (j=startingImage-1; j<max; j+=increment){

      		   open(dir+list[j]);


//selects the frames that were traced in the particle tracking. all other frames that were not tracked return=0 
function indexOfArray(array, value) {
    count=0;
    for (a=0; a<lengthOf(array); a++) {
        if (array[a]==value) {
            count=a;
        }
    }
        return count;
}  // ends function

        in=indexOfArray(cref,j+1);

        //print(j, in, inpast);
           
           if (in > 0){   // selects for the tracked frames
      
                widthx=  -((xref[(in)])-xorigin);
                heighty= -((yref[(in)])-yorigin);       
                inpast=in;   	 //  this index will be used to translate frames that were not tracked. 
               // }
                print(j+1, in, cref[in], xref[in],yref[in]);
               
           } else{      // if they were not tracked it uses the last tracked frame and translates as the last frame was translated.
           	     widthx=  -((xref[inpast])-xorigin);
                 heighty= -((yref[inpast])-yorigin);     
                 print(j+1, inpast, cref[inpast], xref[inpast],yref[inpast]);
               }
           	 
          run("Translate...", "x=widthx y=heighty interpolation=None slice");
          //run("Apply LUT");
          names=File.nameWithoutExtension;
          saveAs("Tiff", destination + names+"-registered");   // saves the individual translated frame
          close(names+".tif");
          close( names+"-registered"+".tif");
          ok = File.delete(destination + names+".tif");

        
          }    //this finish the translations for all movies. 

// now the individual translated frames will be put together in a stack
run("Image Sequence...", "open="+myDir+" number=n starting=1 increment=1 scale=100 file=[] sort");
saveAs("Tiff", path1+"-registered");

//Delete the files and the directory
 list2 = getFileList(myDir);
  for (i=0; i<list2.length; i++)
  ok = File.delete(myDir+list2[i]);
ok = File.delete(myDir);
  if (File.exists(myDir))
      exit("Unable to delete directory");
  else
    //  print("Directory and files successfully deleted");
      print("Adjusted-contrast movie saved as:" );
      print(path1+"-registered");
