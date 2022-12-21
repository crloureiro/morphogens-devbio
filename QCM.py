from __future__ import print_function
from ipywidgets import interact, interactive, fixed, interact_manual, Layout
import ipywidgets as widgets
from IPython.core.display import HTML
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
#import pandas as pd
from math import *


# Define list of possible answers
def choice_list(choices,descript,booleen,val):
    b = widgets.RadioButtons(
                options=[choice for choice in choices],
                value=val,
                description=descript,
                disabled=booleen,
                layout={'width': 'max-content'})
    return b

# Define buttons for possible answers
def button(descript,booleen,icon,style):
    b = widgets.Button(
                description=descript,
                disabled=booleen,
                button_style=style, # 'success', 'info' or 'warning'
                tooltip='description',
                icon=icon)
    return b

# Questions function
def QCM(x,q,a): 
    
    # Display question
    display(HTML(q))
    
    # Define choices
    choices=[x[i] for i in x]
    correct = x['correct']   
    buttons = choice_list(choices,'',False,choices[0])
    
    # Define buttons
    check=widgets.Button(
                  description='Check answers',
                  disabled=False,
                  button_style='', # 'success', 'info' or 'warning'
                  tooltip='description',
                  icon='check')
    
    # Define possible solutions
    check = button('Check answers',False,'check','')
    succes = button('Congratulations!',False,'','success')
    fail = button('Sorry, try again!',False,'','warning')  
    
    # Display choices & buttons
    display(buttons)
    display(check)
    
    # Functions for user interaction
    def on_button_clicked1(b):
        choice = buttons.value
        if choice == correct:
            check.close()
            buttons.close()
            display(choice_list(choices,'',True,choice))
            display(succes)
            display(HTML(a))
        else:
            check.close()
            display(fail)
            
    def on_button_clicked2(b):
        choice = buttons.value
        if choice == correct:
            fail.close()
            buttons.close()
            display(choice_list(choices,'',True,choice))
            display(succes)
            display(HTML(a))
    
    # Check user answers and react accordingly
    check.on_click(on_button_clicked1)
    fail.on_click(on_button_clicked2)

    
# Questions-q & possible answers to each question-c & extra information about the correct answer-a
# Question 1
q1='Question 1: Let´s imagine that Bicoid was targeted during oogenesis using RNAi. What consequences do you think this would have in the oocyte?'

c1={'f1':'Change the position of the nucleus in the cell', 'correct':'Inhibit bicoid by degrading its mRNA', 'f3':'Produce proteins that altered the gradient', 'f4':'Desestabilize the microtubules' }

a1='<FONT COLOR="GREEN">Correct. siRNAs are small RNAs leading to specific knock-off of genes. A siRNA is a small RNA of 21 nucleotides that binds to its targeted mRNA and leads to its mRNA degradation. This results into the absence of the protein of interest. siRNA screenings are widely use in biology to discover new regulators of different genetic pathways. </font>'


# Question 2
q2='Question 2: What feature of the system could be changed to ensure that the gradient diffuses sufficiently fast to reach the posterior part in due time?'

c2={'f1':'Increase translation efficiency', 'correct':'Establish the gradient before translation', 'f3':'Decrease the protein size to improve its diffusion' }

a2='<FONT COLOR="GREEN">Correct. Establishing an mRNA gradient would reduce the time that is necessary to estabslish the protein gradient.  </font>'


# Question 3
q3='Question 3: How could we test this hypothesis? How could we assess bicoid mRNA distribution?'

c3={'correct':'Use a complementary mRNA tagged with a fluorophore', 'f2':'Antibody labeling of targeted mRNA', 'f3':'Integrate a gene encoding a fluorescent protein in the genome' }

a3='<FONT COLOR="GREEN">Correct. The fluorescent probe will bind the bicoid mRNA to form a double strand mRNA which would allow us to see where this mRNA is present inside the oocyte. We cannot use antibodies as they bind to proteins but not to nucleic acids. Integrating a fluorescent probe in the genome could have worked but would still require analyzing the distribution, for instance, of the gfp mRNA. </font>'


# Question 4
q4='Question 4: What do you observe when changing the Diffusion coefficient?'

c4={'f1':'The initial and final concentrations change', 'correct':'The slopes increase when D decreases', 'f3':'The equilibrium is different regarding D'}

a4='<FONT COLOR="GREEN">Correct. Decreasing the Diffusion coefficient leads to an increase of the slopes of the curves. The curves tend to be linear when D decreases. The equilibrium is independent of D, as it only depends on the initial concentration.  </strong>.</font>'


# Question 5
q5='Question 5: What do you think this coefficient $k$ represents?'

c5={'f1':'Cell Division' ,'f2':'Actin dependent transport of proteins' , 'correct':'Protein Degradation' }

a5='<FONT COLOR="GREEN">Correct. $k% is the degradation coefficient. Proteins are transient structures that are degraded constantly within the cell. A change in cell division could affect the gradient, but not the concentration; however, extra cells would result in a larger tissue. An active transport of proteins could be possible; however, actin filaments are not responsible for such a transport.  </strong>.</font>'


# Question 6
q6= 'Question 6: What levels of morphogen do we need to activate a gene promoter for which the transcription factor has a low affinity?'

c6={'f1':'Low morphogen levels only', 'correct':'High morphogen levels only', 'f2':'Both low and high morphogen levels'}

a6='<FONT COLOR="GREEN">Correct. The transcription factor can only bind to the promoter region when it is present at high concentration. Therefore, the effector will be activated only at high concentration of the morphogen.  </strong>.</font>'


# Question 7
q7= 'Question 7: Same question but for high affinity?'

c7={'f1':'Low morphogen levels only', 'f2':'High morphogen levels only', 'correct':'Both low & high morphogen levels'}

a7='<FONT COLOR="GREEN">Correct. Since the binding affinity of the transcription factor is high, it can bind to the promoter of the effector either at high or low concentration. This means that the effector would be activated at low and high concentrations of morphogen.  </strong>.</font>'


# Question 8
q8= 'Question 8: Will cell competence impact the gradient shape?'

c8={'f1':'Yes','correct':'No'}

a8='<FONT COLOR="GREEN">Correct. As the presence or absence of the receptor only impacts how cells perceive the morphogen signal, the differential competence of the cell will not impact the gradient shape. </strong>.</font>'


# Question 9
q9= 'Question 9: Do you expect that the size of the identity domain labeled as C (red domain in previous plot, lacking activation of the effectors) necessarily increases?'

c9={'f1':'Yes','correct':'No'}

a9='<FONT COLOR="GREEN">Correct. The size of the C fate domain will increase only if cells that are exposed to a sufficient morphogen concentration to activate effectors become "non competent".</strong>.</font>'


# Question 10
q10= 'Question 10: Given the same activation threshold for domain, what do you expect to happen with a longer exposure t2>t1?'

c10={'correct':'Domain A will become larger.','f1':'Domain A will remain the same','f2':'Domains A will become smaller'}

a10='<FONT COLOR="GREEN">Correct. If we suppose that all the cells are equally competent and the activation threshold for specific factors do not change, a longer exposure of morphogen will increase the overall signal level which is perceived by the cells. Therefore, the identity domain labeled as A will become wider since a larger part of the tissue is exposed to a morphogen concentration higher than the activation threshold for the A fate. </strong>.</font>'


# Question 11
q11='Question 11: Which of the following parameters does not impact the gradient shape?'

c11={'f1':'Initial concentration', 'f2':'Degradation rate of the morphogen', 'correct':'Cell competence', 'f3':'Diffusion coefficient'}

a11= '<FONT COLOR="GREEN">Correct. The cell competence has an impact on how cells respond to morphogen signaling, not on the gradient itself.</strong>.</font>'


# Question 12
q12='Question 12: What are the parameters that influence the fate cells with adopt?'

c12={'f1':'Distance source', 'f2':'Cell competence','f3':'Exposure duration', 'correct':'All of above'}

a12='<FONT COLOR="GREEN">Correct. By definition, morphogens act differentially on cells depending on their concentration. However, it became clear that differential interpretation by the cells of morphogens depends not only on the <strong>level of signaling</strong> i.e. the distance from the source but also on the signaling dynamics, in particular the <strong>duration of the signal</strong>. Finally, the presence of the morphogen is not enough to induce a response and require the <strong>specific competence of receiving cells</strong>.</font>'