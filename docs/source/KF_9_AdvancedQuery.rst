=============================
Advanced Querying
=============================


ANDs and ORs
**********************************

Once you've built a multiple filter query, your query box should look something
like this:

.. figure:: ./images/KidsFirstPortal_20.png
   :align: center

   **Multiple Filters**

Note that all of your filters are automatically connected by "AND"s, but that
might not always be what you want. If I wanted to see the participants that meet
at least one of my filters, but not all of them, I can click on an "AND" and change
it to "OR", which dramatically changes my results:


|pic21| AND vs OR |pic22|

.. |pic21| image:: ./images/KidsFirstPortal_21.png
   :width: 45%

.. |pic22| image:: ./images/KidsFirstPortal_22.png
   :width: 45%

Note that when you change one AND/OR selection, it changes *all* of the AND/OR
boxes in that query line.

Joint Queries
**********************************

What if our research question was whether children diagnosed with the same brain
tumors have significantly different clinical presentations based on their age
at diagnosis, and we want to compare very young children to teens.

So, we want to add a second filter for "Age at diagnosis" for 15-20 year olds.
However, we can't click on that bar in the plots, and we also don't have the
option to add a second version of that filter in the Filters, we can only change it.

To build a cohort like this, we need to make a more complicated overall query,
where we can have an OR statement about age, but AND statements for our other filters.


have both an 'OR' statement for age and
the AND statem

.. admonition:: Challenge:
   :class: exercise





    * The most expedient way to do this is to click the **Quick Filters** button. Let's try that now.
    * Choose some different options from the **Quick Filters** category and watch how your **Cohort Results** change.
    * To further refine your **Cohort Results** select more categories along the navigation bar.
    * Clicking on the graphics also adds filters to the data
    * Need to cover how to remove the filter
