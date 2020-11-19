# Setting up Terra account


## Linking GCP billing account to Terra billing project

individual account:
- Following the creation of the billing account, they add Terra (terra-billing@firecloud.org) as a “Billing Account User” to link the accounts and complete setup. However, we expect that most users of Common Fund data are not individuals who are personally financing their research, and would instead need to use the more onerous G Suite organization set up.

central billing:
- First, Jon would create an Organization in G Suite, which would make him the “Resource Manager Organization Administrator”. Once the organization has been created, it will automatically have the additional roles of “Resource Manager Project Creator” and “Billing Account Creator”. Jon can then add Janice as either a User or Admin to the G Suite organization. Jon must be an Admin on the G Suite organization to access account management pages on the GCP console (for e.g. to change settings, add/change member permissions, add members). However, since the administrator has the credit card information, Janice has to set up the billing account. Janice would then go to the GCP console to create the Google Billing Account and link a credit card to it. Janice must add Jon to the Google Billing Account as a “Billing Account Administrator” to give him the ability to edit the billing account and add members/permissions to the billing account, including the ability to add Terra as a “Billing Account User”.

## Basics of Terra workspaces

- don't know how much detail, probably cover the basics for tutorial of uploading data and the workflows page
- otherwise link to the support center docs


- set up a workspace, upload data here. do workflow stuff in next page (custom-workflow-terra.md)
