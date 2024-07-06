# Welcome!!

# Thank you for checking my page

My research is on using dynamic Voronoi diagrams to optimize offensive efficiency of women’s lacrosse. In order to achieve this, we track the position of players as a function of time using Python’s Open CV and packages like YOLO and Media Pipe. The data we used was collected from Johns Hopkins Women’s Lacrosse. 
First and foremost, a geometric transformation of the field from side view to bird-eye view has to be determined. The process involves detection of the edges from the video using pixels coordinates and the actual size of the field to create a transformation. Second step is to detect the players using YOLO and Media Pipe and applying the already defined geometric transformation to change the pixel coordinates of the players into actual (x, y) position as a function of time. Lastly, the actual position of the players is used to create Voronoi diagrams which can tell us more about the patterns that emerge right before a successful offensive event like goal or high-quality shot.
We weren’t capable of attaining the real (x, y) coordinates of the players because of the low quality of the video that didn’t enable us to run object detection efficiently. Therefore, we never reached the step of making dynamic Voronoi diagrams. As a result, we started looking into other different ways of tracking lacrosse players on the field by using sensors like wearables.
