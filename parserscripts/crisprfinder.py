"""
Current non-native dependencies in selenium package. If desired to have
script run 'headless' (without any browser opening), then need
additional python selenium PhantomJS dependency.

Instructions for headless operation: 
    - replace 'driver = webdriver.Firefox()' line with
      'driver = webdriver.PhantomJS()'

Note: accessing the CRISPRfinder tool online for many files back to
back may block ip address if website is setup to prevent automated
webcrawling. A potential work-around is to add a 'wait' or 'sleep'
function in between individual calls (i.e. for 50ms).
"""

import os
import sys

from selenium import webdriver

genome_input = sys.argv[1]  # .fasta file format required
crispr_results_output = sys.argv[2]  # textfile with CRISPR results
# path of input/output files is relative to python script location

driver = webdriver.Firefox()
driver.get("http://crispr.i2bc.paris-saclay.fr/Server/")
assert "CRISPR" in driver.title
driver.find_element_by_name("fname").send_keys(
    os.getcwd() + '/' + genome_input
)
driver.find_element_by_name("submit").click()

crispr_results = driver.find_element_by_xpath("//div[@class='content']").text

with open(crispr_results_output, 'wt') as f_handle:
    f_handle.write('\n'.join(crispr_results.split("\n")[3:]))

driver.close()
