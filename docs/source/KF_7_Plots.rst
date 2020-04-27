=============================
Filter Buttons
=============================


.. figure:: ./images/KidsFirstPortal_14.png
   :align: center

   **Explore Data Filters**

To use the filter buttons, click on one, choose from the options in that button,
and click `apply`.

.. admonition:: Exercise:
    :class: exercise

    Use the `Clinical` button to filter the data to include only samples where
    the patient was between 0 and 5 years old when they were diagnosed.

.. hidden-code-block:: python
    :starthidden: True
    :label: --- SHOW SOLUTION TO EXERCISE ---

    Click |H| `Clinical` then `Age at Diagnosis`. Choose `between` in the first
    dropdown, check `year` in the middle, and then type the numbers 0 in the `from`
    box and 5 in the `to` box. Click `Apply`.

If you applied this filter successfully, you should see something like the following:

.. figure:: ./images/KidsFirstPortal_15.png
   :align: center
   :figwidth: 60 %

   **Age at Diagnosis between birth and 5 years old**


Several things about the page have now changed. First, there is now information
in our queries box

.. figure:: ./images/KidsFirstPortal_16.png
   :align: center

   **Queries Box**

It now says what query we are currently looking at, and how many participants are
in our query. Note that it says we chose ages between 0 and 1826.25, which
means it automatically calculated 0-5 years in days. If we wanted to change that
number, we could do so in this bar by clicking the small arrow next to the age range.
A pop up with the same selection box we saw before will appear, and let us change
our selections:

.. figure:: ./images/KidsFirstPortal_17.png
   :align: center
   :figwidth: 50 %

   **Editing Queries**

Cohort Results plots
**********************************

The other big change on this page is that all of the plots have dynamically
updated, and this section has been renamed from "All data" to "Cohort Results
for Query":


.. figure:: ./images/KidsFirstPortal_18.png
   :align: center
   :figwidth: 80 %

   **Cohort Results Plots**

The most notable difference here, is that the distribution in the "Age at Diagnosis"
plot in the center bottom has narrowed dramatically, which is expected given our
query. All of the plots shown in this section respond dynamically to the query
buttons, but they also can be used to build and refine queries, or to see extra
information. For instance, if you hover your mouse pointer over a bar in the
'Experimental Strategy' plot of the "Available Data" section, you'll get the count
of that bar and its x-axis label:


.. figure:: ./images/KidsFirstPortal_19.png
   :align: center
   :figwidth: 60 %

   **Extra Plot Information**

You can continue to filter and refine your cohort by clicking on the bar in a
chart or slices in a pie chart. Note that each time you make a selection in a plot,
that filter is added to your query box above.


.. admonition:: Note:
   :class: tip

   When we ran this query in early 2020, Kids First had data on 222 Participants who
   were diagnosed between the ages of 0 and 5. If you are doing this tutorial later,
   you may see more or different data.


.. admonition:: Exercise:
    :class: exercise

    Using the **Quick Filters** and/or plots further filter your query to only
    include data from the Pediatric Brain Tumor Atlas, and non-Hispanic or Latino
    patients
