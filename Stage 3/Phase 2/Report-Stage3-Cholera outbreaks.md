<!--StartFragment-->

## **Cholera outbreaks Data Analysis Report**

Authors (@slack): Merna Raafat Salem (@MernaSalem28), Sarani Nativel-Santerne (@Sarani)



## **Project files**

### **Summary:**
<https://github.com/MernaSalem/hackbio-cancer-internship-/blob/main/Stage-3/Summary-Stage3-Cholera%20outbreaks.md>

### **R script:** 
<https://github.com/MernaSalem/hackbio-cancer-internship-/blob/main/Stage-3/cholera%20outbreaks%20-Stage3.R>

### **Visualizations:**
<https://github.com/MernaSalem/hackbio-cancer-internship-/tree/main/Stage-3/Visualizations>

Note : for interactive visualizations run the code and they will be downloaded to your pc


## This report outlines the outbreaks and key insights extracted from the WHO Cholera Outbreak Data since 1949.

Cholera is a severe diarrheal disease caused by the bacterium _Vibrio cholerae_, leading to significant morbidity and mortality, especially in regions with inadequate sanitation. Since its first recorded outbreaks, cholera has emerged as a global health threat, prompting ongoing surveillance and response efforts.


## **Data Preparation and Cleaning**

1. Data: The WHO Cholera outbreak Data was loaded into Rstudio

2. Unique and Most Repeated Values: they were identified to understand the structure of the dataset and assess the distribution of key variables.

3. Checking Missing Values and duplicated data: A detailed inspection of missing values across different columns but turned out there was no duplication or missing values.

4. data types are adjusted to ensure numerical consistency 

    Converting **deaths** from **chr** to **integers**, 

    Converting **fatality** rates from **chr** to **double** .


## **Interactive Charts for Enhanced Understanding**

Within the code, you’ll find all charts are interactive to provide deeper insights and enhance your understanding of the Cholera outbreak since 1949,

Ex: fatality rate, deaths and cases over time( line plots ) , number of cases and deaths for the countries …ect. 





## **Key Insights and Summary :** 

**1- Cases caused by cholera :**


![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXduuSmBsY4CzdXTHqQKdhp8XV3K_lVF_Lq4CNWW5HxF5As3ydHMhtJvqIZs4V2JolUMgZJf3s6ZiY6GKHWvon9m1xVFmo6LwsOjmKjIqCN8ka_DgUvysbN0eIourhQ6aFNg5o2PkOVbtxqPAw2Bvc4ficE9?key=9irNdePrIi3O8acveYP2sA)



**1.1 Top 3 countries in recorded number cases:**

- India 1.36325 M

- Haiti 795.794 K

- Peru 736.195 K


![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXed4gV33Ppt-tUqK3QGvMBXEDLKs5EP2QzistJiXQUXudjGZtu28hI4UlcHxo71um1dtx4enq8XoLamExJIhI0G9C8_1lQFZbOPreVwkEb0FYk5CqJcO38GotlHhUpsiQbgN9GaYjmEPF-jRoG-SytondH8?key=9irNdePrIi3O8acveYP2sA)


**1.2 Countries that recorded the lowest number of cases:**

- Lesotho, Mauritius, Slovenia recorded zero cases .

- Hungary, Ireland, Samoa, Uzbekistan recorded only one case .

- Bahamas, Estonia, Greece recorded two cases .

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcExK9YObjI8xMbTxhkB6PzEau_84XAfP6x_q0nB8Jj-JIVcjysEV1g95oHaGnSrxf0abUCrZSQuPJn2zFe7H5WzUCyzFAdUUCUMNelzKfCMCO00mCi3cCPV8KwXKBy_Zq4saJsoZ6Qni9Am1xT9xEse-NL?key=9irNdePrIi3O8acveYP2sA)

**1.3 The country with highest number of cholera Cases**

- In 2011 Haiti  had the highest  recorded number of cholera Cases with 340311 Cases .

**1.4 First Country to Have Highest number Cases**

- Thailand was the first country to record the highest recorded cases as in 1949  the cases were 16 cases .

**2- Deaths from cholera :**

\
![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcOqFiw0GfhEonhTAAMS7QrgB5f81FVVoOsMAZoI23q1q_E_yqmi5_RxQ5IWFKE4hv3A4PVYP6t3K-QM1ggua4OQJzJlmkriZj2b-fDT_ZFlb999nipGSqm9kwxg6iRnrBIWqRWdnE2Z3ofVT5U53RQcbCS?key=9irNdePrIi3O8acveYP2sA)



**2.1 Top 3 countries in recorded number deaths:**

- India 509.438 K

- Bangladesh 137.429 K

- Indonesia 30.6 K

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXe5MbCFNuGZPMAqjs-buhpOoV66B0kwoGOGHEb-3_J7jvPpddRxx7BaMNCwSAPnqHEb9Y9Rt_OwhUYhY4kQO4umd2UZ-DRup9ii_XyRReEuw-lwqrigLiH7MaVtzRR8ismvp3_td0RJvwrAMa3yY4jFCaWQ?key=9irNdePrIi3O8acveYP2sA)

**2.2 Country with highest number of cholera deaths**

- In 1953 India had the highest  recorded  number of deaths which was 124227 deaths .


**2.3 First Country to Have Highest number deaths**

- Thailand was also  the first country to record the highest recorded deaths which was  2 deaths


**3- Fatality rate** ![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcDiH-ZqDVAOKvUDVqSH8Ojz5zuD8Oy994F1bJl6hiGusLMF0-hYdtVTEkAiO2w44-IV4qndzO-CYFAphy7tUFWJSBt3_4vq25NprgEgaJj003XIe_0YZn31Guh35CbHMIuGFso87U-cCeeX_XR_q0JDOY?key=9irNdePrIi3O8acveYP2sA)



**3.1 Top 3 Countries with the Highest fatality rates  over time**

- Bangladesh 1201.33

- India 969.51

- Myanmar 788.98

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXf17vLlQGDfEeIqaGgFWgeuaLK1gAfv6B4bi5ZWOdSFK85qUbIvTkti7mJlf6dnIA83I96A4exAmyjU9BgHN9j_-CpAm48qrj7jAeLyutmuhvbN1vBRZs4dG_uAgQXwxTDlWi1ZAzsuHSb4HRQzVX0WL6w?key=9irNdePrIi3O8acveYP2sA)

**3.2 Highest fatality rate recorded**

- Highest fatality rate was recorded in Italy in 1998 by  450 % .


**3.3 First Country to Have the Highest fatality rate**

- And the first was in China in 1949 by 100% .

**Below is a visualization showing 40 countries with lowest fatality rate** 

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXc74rjTAatdEYbEQkVTb5SMmFOLe-bHpVMhk2TsUSy2Z3FJQPqYkqHjvzKRMcPfxIKXENaT07kVZxfR5qXD5pIEYnUWCCb_IMlYUTK6QBpIH9SiUOwovnYtoqsZw0yFnPi-rJ4G7KmQVdnuMoCHv5gl_LDx?key=9irNdePrIi3O8acveYP2sA)


## **Conclusion:**

This report highlights key insights from the WHO Cholera Outbreak Data, revealing critical trends in case numbers and fatality rates across various countries since 1949. The data emphasizes the ongoing challenge of cholera, particularly in regions like India and Haiti, and underscores the importance of effective public health interventions to mitigate future outbreaks. Our interactive visualizations enhance understanding of these trends, aiding in better preparedness and response strategies.
<!--EndFragment-->
