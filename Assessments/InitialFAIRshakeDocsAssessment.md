# A review of FAIRshake

FAIRshake is a web-based tool that can be used to assess the FAIRness (Find-ability, Accessibility, Interoperability and Reusability) of other web-based biomedical repositories, tools and datasets collectively called digital objects. FAIRshake offers a visual representation of the digital object’s assessment in the form of an insignia to quickly communicate areas of strengths and weaknesses. In addition, FAIRshake also offers a Google Chrome browser extension and a bookmarklet to easily perform FAIRness assessments of websites. 

We had two postdocs look through the FAIRshake documentation to try to answer questions that we have heard from the DCCs about FAIRshake. For this assessment, they focused on an already completed assessment for GTEx. Here is a list of documents they used:
* The FAIRshake website itself: https://fairshake.cloud/
* Documentation page: https://fairshake.cloud/documentation/
* Youtube tutorials on how to use FAIRshake: https://www.youtube.com/playlist?list=PLdix4HBQmtjxCeDo3Px5GxRgXu4GxLotp
* FAIRshake GitHub repo: https://github.com/MaayanLab/FAIRshake
* FAIRshake "report a bug" page: https://github.com/MaayanLab/FAIRshake/issues
* The FAIRshake cookbook: https://github.com/nih-cfde/the-fair-cookbook/blob/c1da35789e7aff4a89ec27fa41c05cc222c26738/content/recipes/04/fairshake.md
* This nature paper: https://www.nature.com/articles/sdata2018118 (We could not actually access this paper)

Overall, they liked the idea of a tool to access fairness. However, they found several aspects of the system unintuitive or difficult to find information on. Below are several suggestions of ways to improve the user experience so that DCC staff can quickly and easily understand their score(s).


## Motivation and documentation
From casually browsing the website (https://fairshake.cloud/project/11/), it is not clear what FAIRshake is, who should care and how to use it. Most of this information is found in the documentation, but it could be useful to have it on the homepage of fairshake. Currently, the documentation tab is located at the bottom of each page in small font. Both reviewers agree that this button needs to be made visible. We recommend: 
* A tab on the FAIRshake website called DOCUMENTATION and 
* Make the homepage more descriptive: all the information provided in the youtube video tutorial description box can go here. This section should also include the motivation for making an object FAIR. It would be good to know why one should bother spending extra time to comply with FAIR standards.

## Insignia
Each page has an insignia that provides a quick and visual representation of the digital object’s FAIRness score. This is very useful. However, at first glance it is hard to interpret. 
* It needs to be clear that one needs to hover on each square to view the specific score. Please use text to indicate hovering options to the user (e.g.: hovering over individual colored boxes on the insignia reveals the digital object's score on said metric). 
  * Please provide a section called "Insignia and how to use/interpret it" in the DOCUMENTATION tab and link it to every page that sports an insignia. 
  * Please indicate which of the FAIR principles that metric relates to - this can be either within the pop-up box of the insignia itself or in the DOCUMENTATION. For example, one metric is phrased "information describing how to cite the tool is provided". Does this relate to FINDABILITY, ACCESSIBILITY, INTEROPERABILITY, or REUSABILITY? For some metrics, this is a straightforwards answer, but that's not always the case.
* Currently, it is not immediately obvious what happens when we click on the insignia. It takes the user to another page with another insignia and an assessments tab. The insignia could be either linked to the “Insignia and how to use/interpret it” page mentioned above or directly to the analytics page with the graphs. Here our assumption is that the data for the squares are also in the analytics graphs. If this assumption is incorrect, please indicate where the numbers that appear on the pop-up box are coming from.
* It is unclear how the many insignias seen in the table in the assessments page are combined to make one insignia for the overall projects page. Assuming that's how the project insignia is created, this problem could be addressed with additional documentation on insignia and score interpretation.
* For some digital objects, the table in the assessments page has [red and blue squares](https://fairshake.cloud/digital_object/566/assessments/), whereas the overall insignia for that tool has [purple and blue squares](https://fairshake.cloud/digital_object/566/). If this confusion is caused by different shades of the same color being used in the insignia vs. the table, please ensure that the colors are exactly the same (for example, the 'yesbut' responses are in green boxes in the assessment table but there are no green squares in the insignia). If not, please provide a detailed key to the color codes used. This key can be located in the "Insignia and how to use/interpret it" section of DOCUMENTATION.
* In an [assessments page](https://fairshake.cloud/metric/15/assessments/) the insignia next to each project heavily implies that the metric is a composite score, when in fact, it is just one of the boxes in the insignia, and the rest of the insignia has nothing to do with this page. This makes it very difficult to interpret the score, as it makes it seem more complicated than it is. We recommend completely removing that column, as it doesn't seem to offer any useful extra information over the final column.

## Analytics page
Similar interpretation problems as the insignia. We recommend a section/tab called "Analytics" in the main DOCUMENTATION tab
* The graphs are not intuitive - there needs to be a key (detailed in the DOCUMENTATION and linked to each analytics page) about what those numbers mean and where they come from. For example: the bars have yes/no responses but the pie chart has only one category (see figure below). Do the numbers on the pie chart correspond to the yes or no responses? In the pie chart below, what does 75% for FAIRshake data rubric mean? Is this the FAIRness score? What is the yesbut category in the bar graphs and how is it calculated?
* Please use text to indicate hovering options to the user (e.g.: hover over graphics with mouse pointer to get more information)

![](https://i.imgur.com/VtFrdoq.png)

* The rubrics tab lists all the metrics used for an assessment. However, there are few details. How many metrics are there? How does it differ for different types of digital objects? We recommend having this information easily accessible to the user, perhaps in the DOCUMENTATION tab. We think it would be useful to know: 
  1) How many metrics our objects are being assessed on 
  2) What metrics are used for each type of digital object 
  3) How FAIRshake finds the information to calculate those metrics
  4) Which metrics correspond to “Findability”, “Accessibility”, etc. If my digital object needs to improve on “Findability”, how do I know which metrics contribute to those scores?

* Additionally it is unclear where one finds the FAIRness score for a project or digital object. Does FAIRshake provide one score per project/digital object/metric? If there is an overall project score, it would be very helpful to show that score on the first page you click on for each project. If there is no overall project score, how do you compare accross projects?  

## Review Criteria
In several instances, we found that different pages of the GTEx website received different ratings, despite being basically the same. For instance in the [Tutorials for the tool are available on the tool’s homepage assessment](https://fairshake.cloud/metric/15/assessments/), 4 out of 5 of the GTEx tools fail, however all of them have documenation featured prominently at the top of the page, and all of them have help and examples on the specific page linked from that assessment. 
As a viewer of this assessment, it's not at all clear what 'counts' as a 'tutorial' here, and so, it's equally unclear what GTEx could do to improve these scores. We suggest providing more, and more descriptive text about what the assessment actually checks, and what the cutoffs are. We also suggest that the assessment include specific feedback about both the Yes's and No's such that GTEx staff could work out what made their one page pass where their others didn't. 

## Reviewer 1 specific suggestions and questions:

It is unclear whether the data used to make the FAIRshake assessments is stored. When you do an assessment, do all the answers get stored as metadata that accompany the FAIRshake scores? It would be helpful for people trying to understand how the FAIRshake requirements were met to see this information.

There are [repeated rows on the GTEx assessments page](https://fairshake.cloud/project/11/assessments/) for 'GTEx Portal Datasets' and 'Gene eQTL Visualizer'. Is this correct? It is affecting the overall FAIRshake score shown for these dataset/tools. 

If one of the primary resources to FAQs is the [journal article published in Cell](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30345-X), it should be made publicly available. Currently, I can only access it if I requested an interlibrary loan. 


## Reviewer 2 specific suggestions and questions:

The utility and difference between chrome extension and FAIRshake bookmarklet are unclear. On their respective tabs, it could be helpful to provide a short explanation of what each of these tools does and why we should care to download them. If I want to do a fairness assessment, do I need to have both? Also, who should download these plugins?








