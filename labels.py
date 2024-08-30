{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "53710d6b",
   "metadata": {},
   "outputs": [],
   "source": [
    "import blabel"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "1a540784",
   "metadata": {},
   "outputs": [],
   "source": [
    "import os \n",
    "import pandas as pd\n",
    "import click\n",
    "import blabel\n",
    "from blabel import LabelWriter"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "14af1388",
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "'/home/alle/jupyter'"
      ]
     },
     "execution_count": 3,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "os.getcwd()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "id": "49f15d4f",
   "metadata": {},
   "outputs": [],
   "source": [
    "mytemplate= open(\"mytemplate.html\", \"w\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "id": "673efaa5",
   "metadata": {},
   "outputs": [],
   "source": [
    "mytemplate.write('''\n",
    "<span class = 'label'>\n",
    "    {{ sample_id }} <br/>\n",
    "</span>\n",
    "''')\n",
    "# Saving the data into the HTML file\n",
    "mytemplate.close()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "870cca0a",
   "metadata": {},
   "outputs": [],
   "source": [
    "import cssutils"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 71,
   "id": "db4ed31f",
   "metadata": {},
   "outputs": [],
   "source": [
    "mystyle1 = open(\"mystyle1.css\", \"w\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 72,
   "id": "224cfecd",
   "metadata": {},
   "outputs": [],
   "source": [
    "mystyle1.write('''\n",
    "@page {\n",
    "    width: 135mm;\n",
    "    height: 17.5mm;\n",
    "    padding: 0.5mm;\n",
    "}\n",
    "\n",
    ".label {\n",
    "    font-family: Times New Roman;\n",
    "    font-weight: bold;\n",
    "    vertical-align: middle;\n",
    "    horizontal-align: left;\n",
    "    display: inline-block;\n",
    "    font-size: 7px;\n",
    "}\n",
    "''')\n",
    "mystyle1.close()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 23,
   "id": "83479d4f",
   "metadata": {},
   "outputs": [],
   "source": [
    "mylabel = pd.read_csv('sppl_labels.csv')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 24,
   "id": "6219c8c1",
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>genotype</th>\n",
       "      <th>sample_id</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>SPxPL</td>\n",
       "      <td>01_0706_01</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>SPxPL</td>\n",
       "      <td>01_0706_02</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>2</th>\n",
       "      <td>SPxPL</td>\n",
       "      <td>01_0706_03</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>3</th>\n",
       "      <td>SPxPL</td>\n",
       "      <td>01_0706_04</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>4</th>\n",
       "      <td>SPxPL</td>\n",
       "      <td>01_0706_05</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "  genotype   sample_id\n",
       "0    SPxPL  01_0706_01\n",
       "1    SPxPL  01_0706_02\n",
       "2    SPxPL  01_0706_03\n",
       "3    SPxPL  01_0706_04\n",
       "4    SPxPL  01_0706_05"
      ]
     },
     "execution_count": 24,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "mylabel.head()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 25,
   "id": "5dadc451",
   "metadata": {},
   "outputs": [],
   "source": [
    "label_writer = LabelWriter('mytemplate.html',\n",
    "                          default_stylesheets=('mystyle1.css',), items_per_page = 2)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 26,
   "id": "1a92247a",
   "metadata": {},
   "outputs": [],
   "source": [
    "records = mylabel.to_dict('records')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 27,
   "id": "bc710777",
   "metadata": {},
   "outputs": [],
   "source": [
    "label_writer.write_labels(records, target = 'sppl_labels.pdf')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "373b0d39",
   "metadata": {},
   "outputs": [],
   "source": [
    "#go to terminal \n",
    "#use function pdfxup (in pdfjam) as follows #working directory contains the pdfxup\n",
    "#if not yet install pdfjam then: sudo apt install textlive-extra-utils\n",
    "\n",
    "\n",
    "#pdfxup -x 1 -y 12 --landscape -fw 0 --tight-frame -o sppl_labels_sheet.pdf sppl_labels.pdf\n",
    "    #1 column, 12 rows #framewidth = 0 (no frame around each label) #tight frame to remove margins\n",
    "    \n",
    "    #open the output in LibreDraw\n",
    "    #set font size to 18\n",
    "    #extend the area so that there is no margin\n",
    "    #the left column is 2cm away from the left margin\n",
    "    #the right column is 17cm away from the left margin\n",
    "    \n",
    "    #save that as odg and (go to some webpage) convert it back to pdf\n",
    "    #ready to print"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.10.6"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
