import os
import sys
import optparse
import pandas as pd
from bs4 import BeautifulSoup
from urllib.parse import urljoin
from requests_html import HTMLSession

"""
Scrape all the known pulsars from different surveys via the pulsar survey scraper (https://pulsar.cgca-hub.org/search)
Mainly based on https://www.thepythoncode.com/article/extracting-and-submitting-web-page-forms-in-python
Input parameters needed: coordinates, search radius, search reference dm and dm tolerance
"""    

URL = "https://pulsar.cgca-hub.org/search"

# initialize an HTTP session
session = HTMLSession()

def get_all_forms(url):
    """Returns all form tags found on a web page's `url` """
    # GET request
    res = session.get(url)
    # for javascript driven website
    # res.html.render()
    soup = BeautifulSoup(res.html.html, "html.parser")
    return soup.find_all("form")

def get_form_details(form):
    """Returns the HTML details of a form,
    including action, method and list of form controls (inputs, etc)"""
    details = {}
    # get the form action (requested URL)
    action = form.attrs.get("action").lower()
    # get the form method (POST, GET, DELETE, etc)
    # if not specified, GET is the default in HTML
    method = form.attrs.get("method", "get").lower()
    # get all form inputs
    inputs = []
    for input_tag in form.find_all("input"):
        # get type of input form control
        input_type = input_tag.attrs.get("type", "text")
        # get name attribute
        input_name = input_tag.attrs.get("name")
        # get the default value of that input tag
        input_value =input_tag.attrs.get("value", "")
        # add everything to that list
        inputs.append({"type": input_type, "name": input_name, "value": input_value})
    # put everything to the resulting dictionary
    details["action"] = action
    details["method"] = method
    details["inputs"] = inputs
    return details


def run_query(opts):
    """ Runs a query on the pulsar scraper with input parameters"""
    # get the first form
    first_form = get_all_forms(URL)[0]
    # extract all form details
    form_details = get_form_details(first_form)
    
    # the data body we want to submit
    data = {}
    for input_tag in form_details["inputs"]:
        if input_tag["type"] == "hidden":
            # if it's hidden, use the default value
            data[input_tag["name"]] = input_tag["value"]
        elif input_tag["type"] != "submit":
            # all others except submit, prompt the user to set it
            if (input_tag['name']) == "coordinates":
               data[input_tag['name']] = opts.coords
            elif (input_tag['name']) == "lb_or_radec":
               data[input_tag['name']] = "y"
            elif (input_tag['name']) == "radius":
               data[input_tag['name']] = opts.radius
            elif (input_tag['name']) == "dm":
               data[input_tag['name']] = opts.dm
            elif (input_tag['name']) == "dmtol":
               data[input_tag['name']] = opts.dm_tol

    # join the url with the action (form request URL)
    url = urljoin(URL, form_details["action"])

    # Send inputs and retrieve data
    res = session.post(url, data=data)
    soup = BeautifulSoup(res.content, "html.parser")
    return soup


def extract_psrs_from_html(soup):     
    """ Gets all relevant pulsar info into a Dataframe """
    # Retrieve all pulsar info from the query
    if 'No pulsars' in soup:
        print ("No pulsars found with the given parameters") 
        sys.exit(0)
    else:       
       # Initialise pandas dataframe
       columns = ['PSR', 'RA (deg)', 'DEC (deg)', 'P (ms)', 'DM (pc cm^-3)', 'Survey', 'Separation (deg)']
       psr_df = pd.DataFrame(columns=columns)       
       all_psr_info = soup.find_all('span')[-2].find("tbody").findAll("tr")
        
       # Generate Pandas Dataframe
       for i,psr in enumerate(all_psr_info):
           psr_df.loc[i] = [psr.findAll("td")[0].text, 
                            psr.findAll("td")[1].text, 
                            psr.findAll("td")[2].text,
                            psr.findAll("td")[3].text,
                            psr.findAll("td")[4].text,
                            psr.findAll("td")[5].text,
                            psr.findAll("td")[7].text]
                                 
       return psr_df 

if __name__== "__main__":

    # Set options
    parser = optparse.OptionParser()
    parser.add_option('--tag', type=str, help='Unique tag for final csv file', dest="tag",default="MeerKAT")
    parser.add_option('--search_coordinates', type=str, help='Reference coordinates', dest="coords",default="12:08 -59:36")
    parser.add_option('--search_radius',type=str,help='Radius in deg. around coordinates to look for psrs',dest="radius",default="1.2") # IB at L-Band
    parser.add_option('--search_dm',type=str,help='Reference DM to search in pc cm^-3',dest='dm',default="")
    parser.add_option('--search_dm_tolerance',type=str,help='search DM tolerance in pc cm^-3',dest="dm_tol",default="10.0")
    parser.add_option('--output_path',type=str,help='Output path to the csv file',dest="output_path",default=os.getcwd())
    opts,args = parser.parse_args() 


    # Run query with input parameters   
    soup = run_query(opts) 
 
    # Get list of pulsars and metadata from the query
    psr_df = extract_psrs_from_html(soup)
    psr_df.to_csv('{}_known_psrs.csv'.format(opts.tag))
    

