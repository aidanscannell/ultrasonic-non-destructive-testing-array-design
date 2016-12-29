# ultrasonic-non-destructive-testing-array-design
Ultrasonic phased array for operation into the human body.

This project focused on the engineering principles beind non destructive testing and not the development of sofware.

This project involved the design of an ultrasonic phased array for operation into the human body. 
The testing medium consisted of 10 points requiring inspection and a back wall 20mm away in the z direction. 
Huygen’s principle was used to generate a calibration array that enabled appropriate parameters to be selected and 
used for the array simulation. 
The time domain signals were then simulated and used in three post-processing algorithms.

<b>Exercise 1</b>: A 2D simulation was produced to predict the beam profile from a linear array operating into water.

<b>Exercise 2</b>: The linear array designed in exercise 1 was optimised for a target point.
It is desirable for the transducer array to focus at the point of interest and converge the acoustic energy at the desired depth.

<b>Exercise 3</b>: A complete set of time-domain signals received from all transmitter-receiver pairs was simulated for 
when the array was used to inspect the sample.
The transducer elements received a 5-cycle toneburst input signal with a centre frequency defined by the transducer setup.

<b>Exercise 4</b>: The simulated time-domain signals from exercise 3 were processed using three imaging algorithms, 
the total focusing method, plane B scan and focused B scan. 
It was found that an optimised array with a higher frequency was more capable of separating the points on the image, 
indicating better resolution at higher frequencies. 
The array was re-optimised following the approach in exercise 2 and the post-processing algorithms were then implemented.
