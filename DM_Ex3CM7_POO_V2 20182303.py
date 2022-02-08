## -*- coding: utf-8 -*-

# Anthony PRAGASSAM 20182303
# M1 GENIOMHE
# Exercice 3 CM7 POO

# Librairie necessaire pour les sequences
from Bio import SeqIO
from Bio.Seq import Seq
import os
import re

# Libraries necessaire pour la gestion de l'Interface Graphique
from tkinter import Tk, messagebox, Text, END, Scale, Frame, HORIZONTAL, Scrollbar, Menu, WORD
from tkinter.filedialog import askopenfilename
import tkinter as tk
from tkinter import filedialog
from tkinter.scrolledtext import ScrolledText
from tkinter import *
try : #Si le module matplotlip n'est pas installé, la fonction graph ne sera pas proposé à l'utilisateur dans le menu affichage
    import matplotlib
    matplotlib.use('TkAgg')
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    import numpy as np
    from matplotlib.figure import Figure
    printgraph = True
except ModuleNotFoundError:
    printgraph = False

Calcul = False
actual_file = False

w = 5
s = 2
var = 1


def new_windows():
####################################################################################
    def open_file():
        global actual_file
        # Fenetre de selection de fichier
        fileselected = askopenfilename(
            filetypes=(("Text File", "*.txt"), ("Fasta file", "*.fsa"), ("Fasta file", "*.fasta")),
            title="Choose a file."
        )
        seqList = list(SeqIO.parse(fileselected, "fasta"))
        if len(seqList) == 0:
            # Le fichier ne contient pas de Fasta proteique!
            messagebox.showerror("Open Error", "Pas de fasta dans le fichier")
        else:
            # Si on a plus d'une sequence dans le fichier
            # on previent l'utilisateur que l'on prends que la 1er
            if len(seqList) > 1:
                messagebox.showwarning("Open Warning",
                                       "Multi fasta detecte! Seul la 1er sequence sera analysee")  # message d'erreur
            # Recuperation de la 1er sequence
            seq_read = seqList[0]
        peptide = Seq(str(seq_read.seq))

        text.delete(1.0, END)

        text.insert(END, peptide)
        actual_file = peptide

    def save_as():
        global actual_file
        actual_file = filedialog.asksaveasfilename()
        open(actual_file + ".txt", "w").write(textrec.get(0., tk.END))

    def save(*useless):
        if not actual_file:
            save_as()
        else:
            open(actual_file, "w").write(textrec.get(0., tk.END))

####################################################################################
    def demarrer():
        #Affichage page Initialisation

        for widget in root.winfo_children():
            widget.pack_forget()
        labelini.pack()
        textContainer.pack()

    def Ftexte():
        # Affichage page des résultats au format texte

        for widget in root.winfo_children():
            widget.pack_forget()
        labeltext.pack()
        textContainer.pack()
        idx = 1.0

        # Code
        #global hydrograph
        #global peptide
        global w
        global s
        #global labelhydro
        textseq.delete(1.0, END)
        textrec.delete(1.0, END)

        text.tag_remove('found', 1.0, END)

        peptide = text.get(1.0, END).upper()
        expression = r"^[KEDRSQNILMVFWCYAGHPT]+$"
        if re.search(expression, peptide): # verification du peptide

            try :
                w = int(Wentry.get())
                s = int(Sentry.get())
            except ValueError :
                w = 5
                s = 2
            hydro = float()
            hydrotemp = float()
            etoile = ""

            for souschaine in peptide:
                if souschaine in "KEDR":
                    hydro = hydro - 1
                elif souschaine in "SQN":
                    hydro = hydro - 0.5
                elif souschaine in "ILMVFW":
                    hydro = hydro + 1
                elif souschaine in "CY":
                    hydro = hydro + 0.5
                elif souschaine in "AGHPT":
                    hydro = hydro + 0

            for i in range(0, len(peptide), 1):
                fenetre = peptide[i:i + w]
                # print(fenetre)
                for souschaine in fenetre:
                    if souschaine in "KEDR":
                        hydrotemp = hydrotemp - 1
                    elif souschaine in "SQN":
                        hydrotemp = hydrotemp - 0.5
                    elif souschaine in "ILMVFW":
                        hydrotemp = hydrotemp + 1
                    elif souschaine in "CY":
                        hydrotemp = hydrotemp + 0.5
                    elif souschaine in "AGHPT":
                        hydrotemp = hydrotemp + 0
                if hydrotemp >= s:
                    etoile = etoile + "*"
                    hydrotemp = 0
                else:
                    etoile = etoile + " "
                    hydrotemp = 0

            # Affichage
            labelhydro = Label(root, text="Score d'hydrophobicité du peptide: " + str(hydro))
            labelhydro.pack()

            labelespace.pack()

            WSContainer.pack()

            textrec.insert(END, ">>\n" + peptide + etoile)
            highlight = varsur.get()
            if highlight == 1:
                sto = ""
                for j in range(0, len(peptide), 60):
                    fenetrepep = peptide[j:j + 60]
                    fenetrestar = etoile[j:j + 60]
                    fenetrepep.rstrip('\n')
                    sto = sto + fenetrepep + fenetrestar + "\n\n"
                textseq.insert(END, sto)
            if highlight == 2:
                sto = ""
                for j in range(0, len(peptide), 60):
                    fenetrepep = peptide[j:j + 60]
                    #fenetrepep.rstrip('\n')
                    sto = sto + fenetrepep +'\n\n\n'
                textseq.insert(END, sto)

                hydro = 0
                for k in range(0, len(peptide), 1):
                    fenetre = peptide[k:k + w]
                    # print(fenetre)
                    for souschaine2 in fenetre:
                        if souschaine2 in "KEDR":
                            hydro = hydro - 1
                        elif souschaine2 in "SQN":
                            hydro = hydro - 0.5
                        elif souschaine2 in "ILMVFW":
                            hydro = hydro + 1
                        elif souschaine2 in "CY":
                            hydro = hydro + 0.5
                        elif souschaine2 in "AGHPT":
                            hydro = hydro + 0
                    if hydro >= s:

                        idx = textseq.search(peptide[k], idx, nocase=1, stopindex=END)

                        lastidx = '% s+% dc' % (idx, len(peptide[k]))

                        textseq.tag_add('found', idx, lastidx)
                        idx = lastidx
                        hydro = 0
                    else:
                        idx = textseq.search(peptide[k], idx, nocase=1, stopindex=END)

                        lastidx = '% s+% dc' % (idx, len(peptide[k]))

                        idx = lastidx
                        hydro = 0

            #textseq.tag_config('found', foreground='red')
            textseq.tag_configure('found', background='yellow', relief='raised')
            textseq.pack()
            SurContainer.pack()
        else:
            messagebox.showerror("Open Error", "La séquence entrée contient des caractères ne provenant pas de séquences peptidique.")
            for widget in root.winfo_children():
                widget.pack_forget()
            labelini.pack()
            textContainer.pack()

    def Fgraph():
        # Affichage page des résultats au format graph

        for widget in root.winfo_children():
            widget.pack_forget()
        labelgraph.pack()
        textContainer.pack()

        global w
        global s

        peptide = text.get(1.0, END).upper()
        expression = r"^[KEDRSQNILMVFWCYAGHPT]+$"
        if re.search(expression, peptide): # verification du peptide

            try:
                w = int(Wentry.get())
                s = int(Sentry.get())
            except ValueError:
                w = 5
                s = 2
            hydro = float()
            hydrotemp = float()
            etoile = ""
            hydrograph = list()

            for souschaine in peptide:
                if souschaine in "KEDR":
                    hydro = hydro - 1
                elif souschaine in "SQN":
                    hydro = hydro - 0.5
                elif souschaine in "ILMVFW":
                    hydro = hydro + 1
                elif souschaine in "CY":
                    hydro = hydro + 0.5
                elif souschaine in "AGHPT":
                    hydro = hydro + 0

            for i in range(0, len(peptide), 1):
                fenetre = peptide[i:i + w]
                # print(fenetre)
                for souschaine in fenetre:
                    if souschaine in "KEDR":
                        hydrotemp = hydrotemp - 1
                    elif souschaine in "SQN":
                        hydrotemp = hydrotemp - 0.5
                    elif souschaine in "ILMVFW":
                        hydrotemp = hydrotemp + 1
                    elif souschaine in "CY":
                        hydrotemp = hydrotemp + 0.5
                    elif souschaine in "AGHPT":
                        hydrotemp = hydrotemp + 0
                if hydrotemp >= s:
                    etoile = etoile + "*"
                    hydrograph.append(hydrotemp)
                    hydrotemp = 0
                else:
                    etoile = etoile + " "
                    hydrograph.append(hydrotemp)
                    hydrotemp = 0

            # Affichage
            labelhydro = Label(root, text="Score d'hydrophobicité du peptide: " + str(hydro))
            labelhydro.pack()

            colgraph = var.get()
            if colgraph == 1:
                colgraph= "-bs"
            if colgraph == 2:
                colgraph="-b"


            lst = list(range(1, len(peptide) + 1))
            y2 = np.array([s])

            fig = Figure(figsize=(30,4.5))
            plt = fig.add_subplot(111)
            plt.axhline(y=0, color='black')
            plt.axhline(y=s, color='gray', linestyle='--', label="Seuil de prédiction")
            plt.plot(lst, hydrograph, colgraph)
            plt.set_xlim([1, len(peptide)])

            plt.fill_between(lst, hydrograph, y2, where=(hydrograph > y2), color='red', alpha=0.8, interpolate=True,
                             label="Partie de la protéine\ntransmembranaire prédite")

            plt.set_title("Profil d’hydrophobicité")
            plt.set_xlabel("Position dans la protéine (aa)")
            plt.set_ylabel("Score d'hydrophobicité simplifiée")
            plt.legend()

            canvas = FigureCanvasTkAgg(fig, master=root)
            canvas.get_tk_widget().pack()
            canvas.draw()

            toolbar = NavigationToolbar2Tk(canvas, root)
            toolbar.update()
            canvas._tkcanvas.pack()

            WSContainer.pack()
            RDContainer.pack()
        else:
            messagebox.showerror("Open Error", "La séquence entrée contient des caractères ne provenant pas de séquences peptidique.")
            for widget in root.winfo_children():
                widget.pack_forget()
            labelini.pack()
            textContainer.pack()

####################################################################################
    # Root Tkinter
    root = tk.Tk()
    root.geometry("900x700")
    root.title("Exercice 3: ")

    # Menu root
    menuBar = tk.Menu(root)
    ##
    fichier = tk.Menu(menuBar)

    fichier.add_command(label="Nouvelle fenêtre ", command=new_windows)
    fichier.add_command(label="Ouvrir un fichier", command=open_file)
    fichier.add_command(label="Enregistrer sous", command=save_as)
    fichier.add_command(label="Enregistrer", command=save)
    fichier.add_separator()
    fichier.add_command(label="Sortir", command=root.quit)
    menuBar.add_cascade(label="Fichier", menu=fichier)

    affichage = tk.Menu(menuBar)
    affichage.add_command(label="Format texte", command=Ftexte)
    if printgraph == TRUE:
        affichage.add_command(label="Format graph", command=Fgraph)
    menuBar.add_cascade(label="Affichage", menu=affichage)

    root.config(menu=menuBar)

    # Création d'un widget Label de la page d'entrée
    label = Label(root, text='Bienvenue',font=(10))
    label.pack()

    # Création d'un widget Label pour ajouter un espace
    labelespace = Label(root, text='')
    labelespace.pack()

    labelexpli = Label(root, text='''
    Ce programme permet de calculer le caractère d'hydrophobicité d'un peptide lu via un menu\nà l'aide de l'échelle d'hydrophobicité simplifiée suivante :
    1.0 : Asp(D), Glu(E), Lys(K), Arg(R)
    -0.5 : Asn(N), Gln(Q), Ser(S)
    +1.0 : Ile(I), Leu(L), Met(M), Val(V), Phe(F), Trp(W)
    +0.5 : Cys(C), Tyr(Y)
    0.0 : Ala(A), Gly(G), His(H), Pro(P), Thr(T)
    
    Il va aussi calculer les scores d'hydrophobicité  d'une fenêtre qui glisse le long de la protéine et l'affiché en fonction d'un seuil. 
    Ce résultat affiché peut être enregistré au format texte.
    
    Enfin, ce programme peut aussi afficher les résultats sous forme de Diagramme d'hydrophobicité.
    
    ''')
    labelexpli.pack()

    # Création d'un widget Button pour démarrer
    button = Button(root, text="Démarrer", command=demarrer,height= 2,width = 10,font=(4))
    button.pack()

    ## Création d'un widget Label de la page d'initialisation
    labelini = Label(root, text='Initialisation',font=( 6))
    #labelini.pack()


## Création d'une frame général
    textContainer = LabelFrame(root, borderwidth=1, relief="sunken", text="Peptide", padx=9)

    ### Création d'une frame avec le bouton pour lancer le positionement des *
    labelbase = Label(textContainer,
                      text='Saisir le peptide ou importer un fichier Fasta puis cliquer sur "Envoyer":')
    buttontext = Button(textContainer, text="Envoyer", command=Ftexte)

    ### Création d'une frame avec le peptide à analyser
    text = Text(textContainer, wrap=NONE, padx=10, pady=10, height=1, width=60, borderwidth=3, relief="sunken")
    textHsb = Scrollbar(textContainer, orient="horizontal", command=text.xview)
    text.configure(xscrollcommand=textHsb.set)

    labelbase.grid(row=0, column=0, sticky="wn")
    buttontext.grid(row=0, column=1, sticky="nsew")
    text.grid(row=1, column=0, columnspan=2, sticky="nsew")
    textHsb.grid(row=2, column=0, columnspan=2, sticky="nsew")

    #textContainer.pack()

    ## Création d'un widget Label de la page de résultat au format texte
    labeltext = Label(root, text="Scores d'hydrophobicité",font=( 6))
    # labeltext.pack()

    ## Saisie taille fenetre et seuil dans une frame
    WSContainer = Frame(root, borderwidth=1)

    Wlabel = Label(WSContainer, text="w, taille de la fenêtre (par défaut = 5):")
    Wentry = Entry(WSContainer, textvariable=w, width=5)
    Wentry.insert(END, w)

    Slabel = Label(WSContainer, text="     s, seuil d'affichage (par défaut = 2):")
    Sentry = Entry(WSContainer, textvariable=s, width=5)
    Sentry.insert(END, s)

    Wlabel.grid(row=0, column=0, sticky="ew")
    Wentry.grid(row=0, column=1, sticky="nsew")
    Slabel.grid(row=0, column=2, sticky="ew")
    Sentry.grid(row=0, column=3, sticky="nsew")

    ## Création d'une zone de texte affichant la séquence + *
    textseq = ScrolledText(root,  padx=10, pady=10,height=11, width=60, bd=1, relief="sunken")

    ## Création d'une zone de texte iddentique à "textseq" pour la fonction enregistrement
    textrec = Text(root,  padx=10, pady=10,height=10, width=60, bd=1, relief="sunken") #

    ## Création d'un widget Label de la page de résultat au format texte
    labelgraph = Label(root, text="Diagramme d'hydrophobicité", font=(6))
    # labelgraph.pack()

    ###Radio Button des options de l'affichage graphique des points
    RDContainer = Frame(root, borderwidth=1)

    var = IntVar()
    var.set(1)
    R1 = Radiobutton(RDContainer, text="Afficher les points (par défaut)", variable=var, value=1,command=Fgraph)
    R2 = Radiobutton(RDContainer, text="Ne pas afficher les points", variable=var, value=2,command=Fgraph)

    R1.grid(row=0, column=0, sticky="nsew")
    R2.grid(row=0, column=1, sticky="nsew")

    ###Radio Button des options de l'affichage graphique des * et surlignage
    SurContainer = Frame(root, borderwidth=1)

    labelsur = Label(SurContainer, text="Légender le caractère d'hydrophobicité via:")

    varsur = IntVar()
    varsur.set(1)
    R3 = Radiobutton(SurContainer, text="Une étoile * (par défaut)     ", variable=varsur, value=1,command=Ftexte)
    R4 = Radiobutton(SurContainer, text="Un surlignage", variable=varsur, value=2,command=Ftexte)

    labelsur.grid(row=0, column=0, columnspan=2, sticky="nsew")
    R3.grid(row=1, column=0, sticky="nsew")
    R4.grid(row=1, column=1, sticky="nsew")

    root.mainloop()



new_windows()

