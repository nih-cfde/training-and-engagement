## Training and Evaluation Plan v2.0 / CFDE / May 2020

**Updates since March - now, with more assessment.**

**Training Plan v2.0 / CFDE / May 2020** </br>

[The training plan for 2020](#the-training-plan-for-2020)
+ [Introduction](#introduction)
+ [The training plan for 2020](#the-training-plan-for-2020)
    + [In-person training vs online training](#in-person-training-vs-online-training)
    + [Online lesson development approach](#online-lesson-development-approach)
    + [Assessment approach](#assessment-approach)
    
[Specific Activities](#specific-activities)
+ [1. Collaborate with DCCs to build program-specific training materials.](#1-collaborate-with-dccs-to-build-program-specific-training-materials)
    + [Tutorial development:](#tutorial-development)
        + [Tutorial Requirements:](#tutorial-requirements)
        + [Video Course Requirements:](#video-course-requirements)
+ [2. Develop general purpose bioinformatics training materials for the cloud.](#2-develop-general-purpose-bioinformatics-training-materials-for-the-cloud)
    + [Tutorial development:](#tutorial-development-1)
        + [Tutorial Requirements:](#tutorial-Requirements-1)
+ [3. Hold discussion sessions to discuss and generalize use cases for data analysis and integration.](#3-hold-discussion-sessions-to-discuss-and-generalize-use-cases-for-data-analysis-and-integration)
+ [4. Discuss and plan opportunities for in-person hackathons for technology development and use case brainstorming.](#4-discuss-and-plan-opportunities-for-in-person-hackathons-for-technology-development-and-use-case-brainstorming)
+ [5. Develop an overarching assessment program.](#5-develop-an-overarching-assessment-program)
    + [The high level goals of our assessment program are to:](#the-high-level-goals-of-our-assessment-program-are-to)
    + [Metrics we will develop and start to measure between now and December:](#metrics-we-will-develop-and-start-to-measure-between-now-and-December)
    
[Hiring](#hiring)
 
## Introduction

This CFDE-CC training plan lays out our long term goals (~3-5 years) as well as our short term plan of action (for the remainder of 2020). It will be followed by another report in October that provides a progress update on 2020 efforts, assessment and evaluation results, and details the next steps for training.

The goals of the CFDE training effort are fourfold. First, we want to work with specific CFDE Data Coordinating Centers to develop and run DCC-specific and targeted cross-DCC-data set training programs to help their users make improved use of their data. Second, we want to provide broad-based training on data analysis in the cloud, to help the entire CFDE user base shift into a more sustainable long term approach. Third, we plan to work with DCCs to execute hackathons in which more advanced users explore “out of the box” data analysis, in order to help DCCs guide their user experience and pipeline development. And fourth, we expect that broad and deep engagement with a wide range of users will help us identify new and exciting use cases for data reuse that can be brought back to the CF DCCs and the CFDE. **Collectively, our training program will train users in basic bioinformatics and cloud computing, help the DCCs lower their support burden, improve user experience, and identify new use cases for data reuse and data integration within and across DCCs.**

All training materials produced by the CFDE-CC will be made available through the central nih-cfde.org web site, under CC0 or CC-BY licenses, which will allow it to be used and remixed by any other stakeholders without limitations. Assessment and iteration on the materials will be carried out by the CFDE-CC’s training team during the pilot period, which we expect to continue through the end of 2020; we may engage with external assessment and evaluation as our efforts expand.

The CFDE-CC’s training component is run by Dr. Titus Brown and Dr. Amanda Charbonneau; as of June 1st, 2020, we have three training postdocs and two staff training coordinators. The training component is closely integrated with the engagement plan, and we expect training to interface with user experience evaluation and iteration across the entire CF and CFDE as well as use case creation and refinement.

In this training plan, we have no specific plans to interface with training efforts outside the Common Fund. However, we are aware of a number of training efforts with similar goals, including Broad’s Terra training program and ANViL’s training focus. Our approaches to sharing materials and running training are designed to allow these other training efforts to make use of our materials, and the underlying technologies and approaches we are using (see below) are entirely compatible.

## The training plan for 2020

#### In-person training vs online training

Our initial plan was to run a series of in-person workshops during 2020. However, our plan is now pivoting to an online strategy, because of the COVID-19 pandemic still sweeping the world. In particular, we expect there to be no in-person meetings before August, and expect there to be substantial barriers to travel after that. We will re-evaluate our plans in October in our updated training plan.

Online training is very different from in-person training. In our experience, in-person training offers a natural focus for many and can support an extended (~4-6 hrs/day) engagement with materials. Moreover, technology problems on the learner’s side can often be fixed by in-person helpers who have direct access to the learner’s computer.  Finally, the intensity of in-person workshops combines well with the higher cost of travel: in the past we have successfully run many in-person workshops, lasting between 2 days and 2 weeks, where either the instructors or the students travel significant distances to attend the workshop.

Online training requires different affordances. Learner attention span in the absence of interpersonal interaction is much shorter. Remote debugging is possible but much less effective than in-person debugging. And both instructors and learners must manage more technology, including teleconferencing software and chat, often on the same screen size as before. These challenges, among others, have limited the effectiveness of online training efforts, including MOOCs; several studies of MOOCs have shown that most learners drop out of MOOCs quickly, and that the main benefits of MOOCs have been to those who already have experience with the material.

In exchange for these downsides, online training offers some opportunities. By using asynchronous delivery of material, different schedules can be accomodated among the learners, and there is much more time for offline experimentation and challenge experiments. Moreover, online training can offer somewhat more scalability and can potentially be offered more cheaply, since it involves no travel or local facilities.

While we believe we can leverage online training effectively, we will need to experiment with formats, try out a variety of technologies, and put more effort into robust material development to offset the challenges of learner debugging. This may delay some of our previously planned training, but should result in robust training materials that meet our original objectives on approximately the same timeline.

#### Online lesson development approach

We need to transition from current draft materials for an in-person workshop to materials that can be delivered online. Our current plan is to start by breaking lessons up into 5-10 minute video chunks, or “vidlets”, that showcase concepts and technical activities. These chunks can be viewed in “flipped classroom” or offline mode, and will be interspersed with opportunities for virtual attendees to seek technical help, explore their own interests, and ask questions in an individual or group setting.

After our initial revamp, we will deliver each lesson within the training team, and then expand to groups outside our team. Each delivery will be a walkthrough of the entire material set with some users, and will result in an iteration to change the materials to reflect discussion during the walkthrough. After 2-3 iterations are delivered to beta users and CF program members, we will set up a formal registration system and encourage adventurous biomedical scientists to attend half-day sessions over a period of a week or two. We have not determined the exact approach we will use, but we expect to combine Zoom teleconferences, livestreaming, and helpdesk sessions.

During the lesson development and delivery period, we will work closely with each partner DCC to make sure our lessons align with their best practices, as well as conveying any technical challenges with user experience back to the DCCs in order to identify potential improvements in DCC portals.

This lesson development approach is slow and cautious, and provides plenty of opportunity to improve the materials in response to lived experience of both instructors and learners. We expect to be able to develop a new lesson (~8-16 hours of total training) approximately every month based on this approach, although we may develop two lessons in tandem; once we have 3-4 lessons developed, we will switch to offering them to a larger audience for a month or two, and then conduct assessments and evaluate our overall approach, as well as next steps for specific lesson development.

#### Assessment approach

Our assessments for this period will be formative and will focus on improving our impact by better understanding the needs of our learners, areas where our materials can be improved, and techniques for better online delivery of our materials. Assessment will primarily consist of within-training check-points, pre- and post-training surveys, and remote interviews with learners and trainers both before and after training. The results from these assessments will be reported in the October training update and will also be used to develop larger-scale instruments that we can use to standardize summative assessment in 2021. We will also work with DCCs to measure continued use by learners as one of our longer-term metrics.

After each scheduled workshop, we will write a summary assessment of the workshop based on any recorded assessment results, and these will serve as one basis for integrated assessment reports and planning in the long term.

## Specific Activities

#### 1. Collaborate with DCCs to build program-specific training materials.

We will work with DCCs to build training materials that help their current and future users make use of their data sets. Our primary goals here are to (a) create and expand materials for users, (b) offer regular trainings in collaboration with the DCCs, (c) provide expanded help documentation for users to lower the DCC support burden, and (d) work with the DCCs over time to refine the user experience and further lower the DCC support burden.

In the near term, Kids First, GTEx, and LINCS have all expressed interest in working with us on specific training opportunities. We have already connected with KF and GTEx, and prior to the COVID-19 pandemic were working on RNAseq tutorials for the KF/Cavatica platform as well as the GTEx/ANViL/Terra platform. With Kids First, we have already conducted one alpha presentation at UC Davis and are communicating with them about our results. For GTEx and LINCs we have specific plans but are waiting for their respective awards.

#### Tutorial development:

We have begun a Kids First WGS and RNAseq tutorial and worked with KF and Cavatica to improve their interface so that it can be used for training. While we continue to work on materials, we are awaiting the chance to iterate with the KF DCC, which is dependent on their hiring processes. We have started work on GTEx materials, however they also are somewhat dependent on input from GTEx. We anticipate ramping this up as GTEx staff come on line in the summer. We plan to deliver at least two sets of workshop materials to the public between now and the end of August, one each on KF WGS/RNAseq and GTEx RNAseq. The exact timeline will depend on how we recruit participants and how quickly our lesson development proceeds. We will be incorporating different pieces (user-led walk-throughs, video lessons, live virtual sessions) and assessing the materials and configurations to deploy the best possible lesson implementations. 

#### Tutorial Requirements:

+ Persistent, user-led walkthrough documents (~2-4 hours of follow-along material each)
+ Accompanying short videos of difficult sections
+ Materials are clinician- focused (KF only)
+ Materials are data scientist- focused (GTEx only)
+ Materials available at nih-cfde.org web site, under CC0 or CC-BY licenses
+ Lessons align with DCC best practices
+ Online community space for learner engagement
    + Formal registration system
    + Code of Conduct
    + Moderator Group
    + DCC and other expert volunteers to answer questions
+ Promotion of materials
+ Promotion of online community
+ Assessment
    + Materials contain breaks for checking understanding/formative assessment.
    + Pre-training surveys on prior knowledge on data sets and techniques, specific learning goals, and self-confidence;
    + Post-training surveys on improved knowledge, learning goals, tutorial format and content, and use case gaps in the training materials.
    + Conduct Remote interviews with learners both before and after training
        + Secure any approvals for human data collection
        + Collect contact information from learners

We have plans with LINCS to help with materials development for a graduate-level MOOC course called “Data Science with LINCS Data” to be offered as an online video series. LINCS has already developed draft materials for this course, but has not had time to test or record them. They would like our team to run through the materials and identify areas that need improvement, and to help them to make those changes. While they would like their own staff to present and record the materials, they have also asked for our assistance in aspects of  post-production such as video editing. Our effort is ready to begin as soon as they are able. 

#### Video Course Requirements:
+ Persistent video lessons
+ Videos are accessible
    + Include written transcripts
    + Include closed-captioning
+ Materials are graduate level, research scientist- focused
+ Materials available at nih-cfde.org web site, under CC0 or CC-BY licenses
+ Lessons align with DCC best practices
+ Promotion of materials
+ Assessment
    + Materials contain breaks for checking understanding/formative assessment.
    + Pre-training surveys on prior knowledge on data sets and techniques, specific learning goals, and self-confidence;
    + Post-training surveys on improved knowledge, learning goals, tutorial format and content, and use case gaps in the training videos.
    + Conduct Remote interviews with learners both before and after training
        + Secure any approvals for human data collection
        + Collect contact information from learners

We will reach out to SPARC to discuss training plans as their funding develops.

#### 2. Develop general purpose bioinformatics training materials for the cloud.

We will develop online training materials for biomedical scientists that want to analyze data in the cloud. Many future NIH Common Fund plans for large scale data analysis rely on analyzing the data on remote-hosted cloud platforms, be they commercial clouds such as Amazon Web Services and Google Compute Platform (GTEx, KF) or on-premise hosting systems like the Pittsburgh Supercomputing Center (HuBMAP). Working in these systems involves several different technologies around data upload, automated workflows, and statistical analysis/visualization on remote platforms.

Since most biomedical scientists have little or no training in these areas, they will need substantial support to take advantage of cloud computing platforms to do large scale data analysis.

#### Tutorial development:

We anticipate providing at least one full set of workshop materials on cloud bioinformatics between now and the end of August. This workshop will consist of a number of different pieces including user-led tutorials, custom made video lessons, and virtual forums and live teaching sessions. The exact timeline and number of workshops will depend on how we recruit participants and how quickly our lesson development proceeds. 

Our “general bioinformatics in the cloud” tutorials are already available for in-person meetings, but need to be updated and revamped to an online format.

For workflows, there are two primary workflow systems in use, WDL and CWL. At least one of these (and sometimes both) is supported by every CF program that uses cloud workflow systems. We will develop initial training materials for data-analysis focused biomedical scientists to make use of these workflow systems, based on our existing workflow materials.

For statistics/visualization, there are two commonly used analysis systems, R/RStudio and Python/Jupyter, that are used by almost all of the CF programs. We already have in-person training material for these systems, and will adapt them to online delivery.

#### Tutorial Requirements:

+ Persistent, user-led walkthrough documents (~4 hours of material each)
+ Accompanying short videos of difficult sections
+ Materials available at nih-cfde.org web site, under CC0 or CC-BY licenses
+ Lessons align with DCC best practices
+ Online community space for learner engagement
    + Formal registration system
    + Code of Conduct
    + Moderator Group
    + DCC and other expert volunteers to answer questions
+ Promotion of materials
+ Promotion of online community
+ Assessment
    + Materials contain breaks for checking understanding/formative assessment.
    + Pre-training surveys on prior knowledge on cloud systems and bioinformatics, specific learning goals, and self-confidence;
    + Post-training surveys on improved knowledge and familiarity with cloud analysis and bioinformatics, improved confidence, learning goals, tutorial format and content, and use case gaps in the training materials
    + Translation of users into help forums, and private communication with volunteers.
    + Conduct Remote interviews with learners both before and after training
        + Secure any approvals for human data collection

#### 3. Hold discussion sessions to discuss and generalize use cases for data analysis and integration.

In tandem with the specific workshops above, we will engage with biomedical scientists who are interested in reusing CF data. This will include members of the CF communities, biomedical scientists who attend our training sessions, and biomedical scientists recruited via social media. These discussions will be used to inform future use case development for data analysis and integration. GTEx in particular is in close contact with their end user community, and has suggested that their user base would be available for engagement.

Assessment efforts for this effort will be generative: as we discover use cases and unmet needs, we will brainstorm ways to assess and evaluate use cases as well as specific metrics to determine if user-identified gaps are being closed.

#### 4. Discuss and plan opportunities for in-person hackathons for technology development and use case brainstorming.

In-person activities such as hackathons and use case brainstorming are incredibly effective ways to develop an understanding of what direction technology needs to move in to enable new data analyses and integrations. While in-person meetings are deferred for the moment, the SPARC and KF DCCs (hosted in Philadelphia) and the Metabolomics DCC (at UCSD) expressed interest in specific events to facilitate technology development. We will connect with all three DCCs to plan events that can be held as travel restrictions ease.

Our assessment modalities for this exercise will focus on identifying individual learning/practice goals in advance, such as “I would like to learn how to analyze collections of SPARC data sets on Amazon Web Services”. We will then evaluate this with a post-event survey asking if “learning goals” were met, combined with targeted interviews with a subset of individual attendees.

#### 5. Develop an overarching assessment program.

#### _The high level goals of our assessment program are to:_

+ Iteratively improve specific trainings as well as overall training program and technology platforms
+ Leave space for open ended conversations to discover new challenges, unmet training needs, negatives and positives about current training effort, and documentation opportunities
+ Iteratively try, adopt, and assess new training technology
+ Provide space for open ended short conversations, to allow us to discover things not covered by surveys

#### _Metrics we will develop and start to measure between now and December:_

+ For online workshops
    + Number of people that show initial interest in training
    + Number of return trainees within a lesson
    + Number of return trainees across lessons
    + Number of trainees that indicate interest in additional as-yet-undeveloped training events
+ For web sites/documentation
    + Site visit metrics
    + Page visit metrics
+ For forums
    + Number of registrations
    + Number of logins
    + Number of posts
    + Number of repeat engagements (e.g. followups to posts)
+ For videos
    + Video watch statistics
    + Video completion statistics
    + Web site hosting stats

We will assess overall confidence metrics for both “is this training potentially relevant/useful based on its description” and self-confidence in actualizing bioinformatics analyses. See https://carpentries.org/assessment/ for some examples.

We will also explore a number of technologies to measure within-lesson engagement and do formative assessment. While asynchronous online training challenges traditional “stop-and-quiz” approaches, low-stakes multiple-choice quizzes can be incorporated into online lessons easily and provide valuable feedback to learners and trainers. Faded examples that learners can fill in on their own time and submit via a common interface can be used to provide feedback asynchronously. More dynamic documentation, supporting both quizzes and executable code, could be used to provide engaging exercises. However, all of these require experimentation and evaluation in order to determine which choices work best within the context of the platforms we choose to host videos and tutorials. This experimentation will be an ongoing part of our training work, and will be reviewed in the October training plan.

## Hiring

As of the beginning of June, we have completed our hiring for training in 2020, and have three dedicated training postdocs and two staff training coordinators.
