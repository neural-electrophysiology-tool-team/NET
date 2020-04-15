# Getting Started with Visual Stimulation
## Data representation
The recorded data can be saved in many formats. In dataRecording.m (located in *NET\\timeSeriesViewer\\dataRecordingObjects\\@dataRecording*) you can find the dataRecording abstract class which defines the interface through which to access the recorded data. It contains 4 important function for reading the data from the files – *getData*, *getAnalogData*, *getTrigger* and *getDigitalData*. Their parameters are documented inside the file.

Now, when we want to import data from a file we need to identify the file’s format and then use the proper class that supports the format. For example, when we want to import data from the “binary“ format, we will use the class *binaryRecording*;

    br = binaryRecording('\\sil2\Data\Turtle\U4_071014\U4_071014_Images3001.bin');

We can then use getData to load the data for specified times and windows, or, if we want to use the timeSeriesViewer GUI to handle that for us and plot the data we can just pass the binary recording to it:

    timeSeriesViewer(br)

## Triggers
Triggers are useful to indicate events that happened through a recording. For example, the time when an image has started appearing to the subject, and the time when the image stopped showing. Using triggers we can associate those events with the cortical data. This is necessary for drawing conclusions based on the data.

To get the triggers from a recording we use *getTrigger*:

    triggers = br.getTrigger()

The triggers are user-defined and change between different simulations. They usually contain lists of start times and end times corresponding to the defined trigger (continuing the example above, a picture has appeared or disappeared) and maybe a trigger for the stimulation, meaning there will be two lists of size 1 that contain the start and end times of the whole stimulation.

After getting the triggers, we can use those triggers to get the cortical data in these times. For example:

    data = br.getData([],triggers{5}(1:10),1500);

Here, in the 5th cell array of the triggers resides the times in which a picture appeared.

## Spike Sorting
After loading the data, we need to identify the spikes in the data. We use a spike sorting algorithm to do so (such as kilosort). That way we can generate a “t-ic” container which holdes information about the spikes:
-   t is an array of the times [ms] which contain spikes
-   ic is a matrix that has 4 rows and N columns (corresponding to N neurons). The rows are: 
    -   channel number (electrode number) 
    -   neuron number 
    -   start index in "t"
    -   end index in "t"

For example – to get the spike times of the first neuron:

    t(ic(3,1):ic(4,1))

### Example of spike sorting
*TODO*

## Raster Plots
A simple visual representation of spikes. It can be generated using *BuildBurstMatrix* (a pre-compiled c function located in *NET\general functions*). For example:

    bm = BuildBurstMatrix(ic, t, triggers{5}(1:10).', 1000);

Where the parameters are as follows:
1.  The “ic” matrix we got from the spike sorting. 
2.  The “t” array of spike times [ms]. 
3.  The times [ms] from which to take a burst window. The amount of times here is corresponding to the number of bursts (or trials) inside the output matrix.
4.  The width of the burst window [ms].
    
The output matrix is of size [trials x neurons x width]. We can then plot the raster plot using *imagesc*, taking one trial at a time:

    imagesc(squeeze(bm(1,:,:)))

