=============================
Exploring Data in the Portal
=============================

Now that all of our accounts are interconnected, let's actually explore the data.


Step 1 Get logged in to the Kids First DRC Portal
**************************************************

If you aren't still logged in to the Kids First Portal, navigate to
`https://portal.kidsfirstdrc.org/ <https://portal.kidsfirstdrc.org/>`_ and use the
credentials you set up previously in this lesson to log in. When you first log in,
you will see this screen:

.. figure:: ./images/KidsFirstPortal_11.png
   :align: center

   **Kids First Dashboard.**


Step 2 Navigate to the Explore Data tab
**************************************************

Click on `Explore Data`

.. figure:: ./images/KidsFirstPortal_12.png
   :align: center

   **Go To Explore Data**

You should now be on a page that looks like this:

.. figure:: ./images/KidsFirstPortal_13.png
   :align: center

   **Explore Data**


***************************************************
Exploring the Data
***************************************************


It is important to get a sense of what data exists before we start filtering
down and designing an experiment. This page gives us a set of six interactive plots
that display the overall data. You can impose filters on the data in two different
ways:

* By using the filter buttons
* By clicking a graph component of any plot

Filter buttons
**********************************

.. figure:: ./images/KidsFirstPortal_14.png
   :align: center

   **Explore Data Filters**

To use the filter buttons, click on one, choose from the options in that button,
and click `apply`. Watch how this changes the plots below, and note that the
title for the plot area changes from "All data" to "Cohort Results
for Query"

.. hidden-code-block:: python
    :starthidden: False

    a = 10
    b = a + 5

Using ``label`` change the toggle text and ``linenos``
to include line numbers:

.. hidden-code-block:: python
    :linenos:
    :label: --- SHOW ANSWER ---

    x = 10
    y = x + 5


..  Exercise:
    Use the `Clinical` button to filter the data to include only samples where
    "Age at Diagnosis" is between 5 and 10 years.


..  Exercise:
    Choose some different options from the **Quick Filters** category and watch
    how your plots change.


    * The most expedient way to do this is to click the **Quick Filters** button. Let's try that now.
    * Choose some different options from the **Quick Filters** category and watch how your **Cohort Results** change.
    * To further refine your **Cohort Results** select more categories along the navigation bar.
    * Clicking on the graphics also adds filters to the data
    * Need to cover how to remove the filter
