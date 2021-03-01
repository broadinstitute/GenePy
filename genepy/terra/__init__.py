# Jeremie Kalfon
# jkobject@gmail.com
# for BroadInstitute CDS
# in 2019

import time
import pandas as pd
import dalmatian as dm
import numpy as np
import os
import re
import signal
from genepy.utils import helper as h
from genepy.google import gcp
import pdb
import subprocess
from gsheets import Sheets


def createManySubmissions(workspace, workflow, references, entity=None, expression=None, use_callcache=True):
  """
  wrapper to create many submissions for a workflow

  Args:
  ----
    workspace: str namespace/workspace from url typically
    references: list(str) a list of name of the row in this entity
    entity: str terra csv type (sample_id...)
    expresson: str to use if want to compute on the direct value of the entity or on values of values
                e.g. this.samples
    use_callcache: Bool to false if want to recompute everything even if same input

  Returns:
  ------
    submission_ids list(str) the submission ids
  """
  wm = dm.WorkspaceManager(workspace)
  submission_ids = []
  for ref in references:
    submission_ids += [wm.create_submission(workflow, ref, entity, expression, use_callcache)]
  return submission_ids


def waitForSubmission(workspace, submissions, raise_errors=True):
  """
  wrapper to create many submissions for a workflow

  Args:
  -----
    workspace: str namespace/workspace from url typically
    submissions: list[str] of submission ids
    raise_errors: bool to true if errors should stop your code

  Returns:
  -------
    list of ids of failed submissions
  """
  failed_submission = []
  timing = 0
  wm = dm.WorkspaceManager(workspace).disable_hound()
  assert submissions is not None
  if type(submissions) is type(""):
    submissions = [submissions]
  for scount, submission_id in enumerate(submissions):
    finished = False
    while not finished:
      done = 0
      failed = 0
      finished = True
      submission = wm.get_submission(submission_id)["workflows"]
      for _, i in enumerate(submission):
        if i['status'] not in {'Done', 'Aborted', 'Failed', 'Succeeded'}:
          finished = False
        elif i["status"] in {'Failed', 'Aborted'}:
          failed += 1
          if i["workflowEntity"]["entityName"] not in failed_submission:
            print(i["workflowEntity"]["entityName"])
            failed_submission.append(i["workflowEntity"]["entityName"])
        elif i["status"] in {'Done', 'Succeeded'}:
          done += 1
      if not finished:
        time.sleep(40)
        print("status is: Done for " + str(done) + " jobs in submission " + str(scount) + ". " + str(timing) + ",5 mn elapsed.", end="\r")
        timing += 1
        time.sleep(20)
        print("status is: Failed for " + str(failed) + " jobs in submission " + str(scount) + ". " + str(timing) + " mn elapsed.", end="\r")
      else:
        print(str(done / (done + failed)) + " of jobs Succeeded in submission " + str(scount) + ".")
  if len(failed_submission) > 0 and raise_errors:
    raise RuntimeError(str(len(failed_submission)) + " failed submission")
  return failed_submission
  # print and return well formated data


def removeSamples(workspace, samples):
  """
  removes a set of samples from a workspace (very usefull when we have linked pairs and pairsets)

  Args:
  -----
    workspace: str workspace name
    samples: list of samples
  """
  wm = dm.WorkspaceManager(workspace).disable_hound()
  try:
    wm.delete_sample(samples)
  except:
    print('we had pairs.')
    pairs = wm.get_pairs()
    pairid = pairs[pairs.case_sample.isin(samples)].index.tolist()
    for k, val in wm.get_pair_sets().iterrows():
      wm.update_pair_set(k, set(val.tolist()[0]) - set(pairid))
    wm.delete_pair(pairid)
    wm.delete_sample(samples)


def uploadFromFolder(gcpfolder, prefix, workspace, sep='_', loc=0,
                     fformat="fastq12", newsamples=None, samplesetname=None, source='U',
                     bamcol="bam", baicol="bai", test=True):
  """
  upload samples (virtually: only creates tsv file) from a google bucket to a terra workspace

  A very practical function when you have uploaded a folder of samples (RNA/WES/...) in google storage
  with some naming convention, and want to have them all well listed in Terra for futher processing
  it can create a sample set.
  for a set of files: gs://bucket/path/to/files


  Args:
  -----
    gcpfolder: a gs folder path
    prefix: str the folder path
    workspace: str namespace/workspace from url typically
    sep: str the separator (only takes the first part of the name before the sep character)
    fformat bambai, fastq12, fastqR1R2 given the set of files in the folder (they need to be in this naming format)
            e.g. name.bam name.bai / name1.fastq name2.fastq / name_R1.fastq name_R2.fastq
    newsamples: DONT USE
    samplesetname: str all uploaded samples should be part of a sampleset with name..

  Returns:
  --------
    the uploaded dataframe
  """
  wm = dm.WorkspaceManager(workspace)
  print('please be sure you gave access to your terra email account access to this bucket')
  if samplesetname is None:
    samplesetname = 'from:' + gcpfolder + prefix
  files = gcp.list_blobs_with_prefix(gcpfolder, prefix, '/')
  if fformat == "bambai":
    if newsamples is None:
      data = {'sample_id': [], bamcol: [], baicol: []}
      for file in files:
        if file.split('.')[-1] in ["bam", "bai"]:
          name = file.split('/')[-1].split('.')[0].split(sep)[loc]
          if name in data['sample_id']:
            pos = data['sample_id'].index(name)
            if file[-4:] == ".bam":
              data[bamcol].insert(pos, 'gs://' + gcpfolder + '/' + file)
            elif file[-4:] == ".bai":
              data[baicol].insert(pos, 'gs://' + gcpfolder + '/' + file)
          else:
            data['sample_id'].append(name)
            if file[-4:] == ".bam":
              data[bamcol].append('gs://' + gcpfolder + '/' + file)
            elif file[-4:] == ".bai":
              data[baicol].append('gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
        else:
          print("unrecognized file type : " + file)
      df = pd.DataFrame(data)
      df = df.set_index("sample_id")
      print(df)
      df["participant"] = pd.Series(data['sample_id'], index=data['sample_id'])
      wm.upload_samples(df)
      wm.update_sample_set(samplesetname, df.index.values.tolist())
    else:
      # TODO: check if each column exists and can be added, else don't add it
      for i, val in enumerate(newsample["file_path"]):
        if val.split('/')[-1].split('.')[1] != "WholeGenome" or val.split('/')[-2] != "bam":
          newsample = newsample.drop(i)
        elif val.split('/')[1] != 'gs:':
          newsample["file_path"][i] = gcpfolder + newsample["file_path"][i].split('/')[-1]
      newsample = newsample.reset_index(drop=True)
      newsample = newsample.rename(index=str, columns={"sample_name": "sample_id", "subject_name": "participant_id", "file_path": "WGS_bam"})
      currfile = ""
      bai = [''] * int(newsample.shape[0])
      # creating an array of bai and adding it to their coresponding bams
      for i in newsample.index:
        currfile = newsample["WGS_bam"][i]
        if currfile.split('/')[-1].split('.')[-1] == "bai":
          bai[int(newsample[newsample["WGS_bam"] == currfile[:-4]].index.values[0])] = currfile
      newsample["WGS_bam_index"] = pd.Series(bai, index=newsample.index)
      # removing original bai rows
      for i in newsample.index:
        currfile = newsample["WGS_bam"][i]
        if currfile.split('/')[-1].split('.')[-1] == "bai":
          newsample = newsample.drop(i)
      newsample = newsample.reset_index(drop=True)
      newsample["sample_set"] = pd.Series([samplesetname] * int(newsample.shape[0]), index=newsample.index)
      newsample.set_index("sample_id", inplace=True, drop=True)
      newsample = newsample[newsample.columns.tolist()[1:] + [newsample.columns.tolist()[0]]]
      newsample = newsample.loc[~newsample.index.duplicated(keep='first')]
      newsample.to_csv("temp/samples.bambai.tsv", sep="\t")
      wm.upload_samples(newsample)
      wm.update_sample_set(samplesetname, newsample.index)
  if fformat in {"fastq12", "fastqR1R2"}:
    data = {'sample_id': [], 'fastq1': [], 'fastq2': []}
    # print and return well formated data
    print(files)
    for file in files:
      if file[-9:] == ".fastq.gz" or file[-6:] == ".fq.gz":
        name = re.split(sep, file.split('/')[-1].split('.')[0])[loc]
        if name in data['sample_id']:
          pos = data['sample_id'].index(name)
          if fformat == "fastqR1R2":
            if "R1" in file:
              data['fastq1'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            elif "R2" in file:
              data['fastq2'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
          else:
            if file.split('.')[-3][-1] == '1':
              data['fastq1'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            elif file.split('.')[-3][-1] == '2':
              data['fastq2'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq 1/2 error", file)
        else:
          data['sample_id'].append(name)
          if fformat == "fastqR1R2":
            if "R1" in file:
              data['fastq1'].append('gs://' + gcpfolder + '/' + file)
            elif "R2" in file:
              data['fastq2'].append('gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
          else:
            if file.split('.')[-3][-1] == '1':
              data['fastq1'].append('gs://' + gcpfolder + '/' + file)
            elif file.split('.')[-3][-1] == '2':
              data['fastq2'].append('gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
      else:
        print("unrecognized file type : " + file)
    df = pd.DataFrame(data)
    df["Source"] = source
    df["participant"] = data['sample_id']
    df = df.set_index("sample_id")
    if not test:
      wm.upload_samples(df)
      wm.update_sample_set(samplesetname, df.index.values.tolist())
    return df


def updateAllSampleSet(workspace, newsample_setname, Allsample_setname='All_samples'):
  """
  update the previous All Sample sample_set with the new samples that have been added.

  It is especially useful for the aggregate task. Can more generally merge two samplesets together

  Args:
  ----
    workspace: str namespace/workspace from url typically
    newsample_setname: str name of sampleset to add to All_samples
  """
  prevsamples = list(dm.WorkspaceManager(workspace).get_sample_sets().loc[Allsample_setname]['samples'])
  newsamples = list(dm.WorkspaceManager(workspace).get_sample_sets().loc[newsample_setname]['samples'])
  prevsamples.extend(newsamples)
  dm.WorkspaceManager(workspace).update_sample_set(Allsample_setname, list(set(prevsamples)))


def addToSampleSet(workspace, samplesetid, samples):
  """
  add samples to a sample set

  will create new if doesn't already exist, else adds to existing

  Args:
  ----
    workspace: the workspace name
    samplesetid: the sample set name
    samples: a list of samples
  """
  try:
    prevsamples = dm.WorkspaceManager(workspace).get_sample_sets()['samples'][samplesetid]
    samples.extend(prevsamples)
  except KeyError:
    print('The sample set ' + str(samplesetid) + ' did not exist in the workspace. Will be created now...')
  dm.WorkspaceManager(workspace).update_sample_set(samplesetid, list(set(samples)))


def addToPairSet(workspace, pairsetid, pairs):
  """
  add pairs to a pair set

  will create new if doesn't already exist, else adds to existing

  Args:
  ----
    workspace: the workspace name
    pairsetid: the pair set name
    pairs: a list of pairs
  """

  try:
    prevpairs = dm.WorkspaceManager(workspace).get_pair_sets().loc[[pairsetid]].pairs[0]
    pairs.extend(prevpairs)
  except KeyError:
    print('The pair set ' + str(pairsetid) + ' did not exist in the workspace. Will be created now...')
  dm.WorkspaceManager(workspace).update_pair_set(pairsetid, list(set(pairs)))

# Gwen's old version - caught some niche conditions made by get_pair_sets()
# that I think may raise errors in the current version.
# yep. the niche condition has to do with the fact that a list isn't hashable
# def addToPairSet(wm, pairsetid, pairs):
#   pairsets = wm.get_pair_sets()
#   prevpairs = pairsets.loc[[pairsetid]].pairs.tolist() # is this always a list of list? I think so.
#   print(type(prevpairs[0]))
#   if isinstance(prevpairs[0], str) :
#     pairs.extend(prevpairs)
#   elif isinstance(prevpairs[0], list):
#     pairs.extend(prevpairs[0])
#   wm.update_pair_set(pairsetid, list(set(pairs)))


def saveOmicsOutput(workspace, pathto_cnvpng='segmented_copy_ratio_img',
                    pathto_stats='sample_statistics',
                    specific_cohorts=[],
                    speicifc_celllines=[],
                    is_from_pairs=True,
                    pathto_snv='filtered_variants',
                    pathto_seg='cnv_calls',
                    datadir='gs://cclf_results/targeted/kim_sept/',
                    specific_samples=[]):
  """
  *WIP* For a workspace containing all omics workflows (CNV/SNV) (like CCLF's) copies all interesting output to a data bucket

  Args:
  -----
    workspace: the workspace name
    pathto_cnvpng: sample col of the CNV plot results
    pathto_stats: sample col of the bam QC results
    specific_cohorts: if provided, will only look for this specific
    speicifc_celllines: if need to rrun on specific cell lines
    is_from_pairs: if we process on pairs or samples data
    pathto_snv: sample col of the snv files
    pathto_seg: sample col of the segment files
    datadir: gs bucket path where to copy the resulting files
    specific_samples: if provided will only look for these samples

  """
  if specific_cohorts:
    samples = dm.WorkspaceManager(workspace).get_samples()
    samples = samples[samples.index.isin(specificlist)]
  if is_from_pairs:
    pairs = dm.WorkspaceManager(workspace).get_pairs()
    pairs = pairs[pairs['case_sample'].isin(specificlist)]
  for i, val in samples.iterrows():
    os.system('gsutil cp ' + val[pathto_seg] + ' ' + datadir + i + '/')
    os.system('gsutil cp ' + val[pathto_cnvpng] + ' ' + datadir + i + '/')
    os.system('gsutil cp ' + val[pathto_stats] + ' ' + datadir + i + '/')
    if is_from_pairs:
      snvs = pairs[pairs["case_sample"] == i][pathto_snv]
      for snv in snvs:
        if snv is not np.nan:
          os.system('gsutil cp ' + snv + ' ' + datadir + i + '/')
          break
    else:
      os.system('gsutil cp ' + val[pathto_snv] + ' ' + datadir + i + '/')


def changeGSlocation(workspacefrom, newgs, workspaceto=None, prevgslist=[], index_func=None,
                     flag_non_matching=False, onlysamples=[], onlycol=[], entity='samples', droplists=True, keeppath=True, dry_run=True, par=20):
  """
  Function to move data around from one workspace to a bucket or to another workspace.

  can also work on dataframes containing lists of paths

  Args:
  -----
    workspacefrom: the workspace name where the data is
    newgs: the newgs bucket where to copy the data in
    workspaceto: if we should have these new samples and columns added to another workspace instead \
    of just updating the same one (usefull to copy one workspace to another)
    prevgslist: if providded, will only move files that are in the set of google bucket listed here
    index_func: *WIP* unused
    flag_non_matching: if set to true and prevgslist is set to some value, will return a list of samples that were not
    matched to anything in the prevgslist
    onlycol: do this only on a subset of columns in terra workspace
    entity: the entity in the terra workspace on which to do this
    droplists: if set to true remove all columns containing list of paths (list of path are not uploaded well in terra)
    keeppath: if set to true, will keep the full object path and just change the bucket
    dry_run: if set to true will not update anything on Terra but just return the result
    par: on how many processor do the gs copy commands.

  Returns:
  -------
    torename: the pandas.df containing the new paths
    flaglist: the samples that were non matching (if flag_non_matching is set to true)
  """
  flaglist = []
  wmfrom = dm.WorkspaceManager(workspacefrom)
  a = wmfrom.get_entities(entity)
  if len(onlysamples) > 0:
    a = a[a.index.isin(onlysamples)]
  print("using the data from " + workspacefrom + " " + entity + " list")
  if len(a) == 0:
    raise ValueError('no ' + entity)
  if onlycol:
    a = a[onlycol]
  todrop = set()
  torename = {}
  print('this should only contains gs:// paths otherwise precise columns using \"onlycol\"')
  for col in a.columns.tolist():
    val = []
    for k, prev in a[col].iteritems():
      if type(prev) is str:
        new = prev
        if newgs not in new:
          if len(prevgslist) > 0:
            for prevgs in prevgslist:
              new = new.replace(prevgs, newgs)
            if flag_non_matching:
              if new == prev:
                flaglist.append(prev)
          if not keeppath:
            new = newgs + new.split('/')[-1]
          else:
            new = newgs + '/'.join(new.split('/')[3:])
        else:
          print("sample " + str(k) + " was already in the new gs")
        val.append(new)
    # IN CASE WE HAVE A LIST
      if type(prev) is list:
        if droplists:
          todrop.add(k)
          continue
        ind = []
        for prevname in prev:
          newname = prevname
          if newgs not in newname:
            if len(prevgslist) > 0:
              for prevgs in prevgslist:
                new = new.replace(prevgs, newgs)
              if flag_non_matching:
                if new == prev:
                  flaglist.append(prev)
            if not keeppath:
              new = newgs + new.split('/')[-1]
            else:
              new = newgs + '/'.join(new.split('/')[3:])
          else:
            print("sample " + str(k) + " was already in the new gs")
          ind.append(newname)
        val.append(ind)
    torename.update({col: val})
    if not dry_run:
      if keeppath:
        h.parrun(['gsutil mv ' + a.iloc[i][col] + ' ' + v for i, v in enumerate(val)], cores=20)
      else:
        gcp.mvFiles(a[col].tolist(), newgs)
    else:
      if keeppath:
        print(['gsutil mv ' + a.iloc[i][col] + ' ' + v for i, v in enumerate(val)])
      else:
        print("mv "+str(a[col].tolist()) +" "+newgs)
  if workspaceto is None:
    wmto = wmfrom
  else:
    wmto = dm.WorkspaceManager(workspaceto)
  torename = pd.DataFrame(data=torename, index=[i for i in a.index.tolist() if i != 'nan'])
  if not dry_run:
    wmto.disable_hound().update_entity_attributes(entity, torename)
  return torename, flaglist


def renametsvs(workspace, wmto=None, index_func=None):
  """
  ################## WIP ############
  only works for one use case
  """
  data = {}
  wmfrom = dm.WorkspaceManager(workspace)
  try:
    a = wmfrom.get_participants()
    data.update({'participants': a})
  except:
    print('no participants')
  try:
    a = wmfrom.get_samples()
    data.update({'samples': a})
  except:
    print('no samples')
  try:
    a = wmfrom.get_pair_sets()
    data.update({'pair_sets': a})
  except:
    print('no pair_sets')
  try:
    a = wmfrom.get_pairs()
    data.update({'pairs': a})
  except:
    print('no pairs')
  try:
    a = wmfrom.get_sample_sets()
    data.update({'sample_sets': a})
  except:
    print('no sample_sets')
  # currently works only for sample, sample
  for k, entity in data.items():
    ind = []
    for i in entity.index:
      pos = val.find('-SM')
      if pos != -1:
        val = val[pos + 1:]
        pos = val.find('-SM')
        if pos != -1:
          val = val[:9] + val[pos + 1:]
      ind.append(val)
    entity.index = ind
    # for all columns of the tsv
    for k, val in entity.iterrows():
      for i, v in enumerate(val):
        if type(v) is list or type(v) is str:
          ind = []
          for j in v:
            pos = j.find('-SM')
            if pos != -1:
              j = j[pos + 1:]
              pos = j.find('-SM')
              if pos != -1:
                j = j[:9] + j[pos + 1:]
            ind.append(j)
          val[i] = ind
        entity.loc[k] = val
    if wmto is None:
      wmto = wmfrom
    if "participants" in data:
      wmto.upload_participants(data['participants'].index.tolist())
    if "samples" in data:
      wmto.upload_samples(data['samples'])
    if "pairs" in data:
      wmto.upload_entities('pair', data['pairs'])
    if "pair_set" in data:
      pairset = data['pair_set'].drop('pairs', 1)
      wmto.upload_entities('pair_set', pairset)
      for i, val in data['pair_set'].iterrows():
        wmto.update_pair_set(i, val.pairs)
    if "sample_set" in data:
      sampleset = data['sample_set'].drop('samples', 1)
      wmto.upload_entities('sample_set', sampleset)
      for i, val in data['sample_set'].iterrows():
        wmto.update_sample_set(i, val.samples)


def findBackErasedDuplicaBamteFromTerraBucket(workspace, gsfolder, bamcol="WES_bam", baicol="WES_bai"):
  """
  If you have erased bam files in gcp with bai files still present and the bam files are stored elsewhere
  and their location is in a terra workspace.

  Will find them back by matching bai sizes and copy them back to their original locations

  Args:
  ----
    workspace: str namespace/workspace from url typically
    gsfolder: str the gsfolder where the bam files are
    bamcol: str colname of the bam
    baicol: str colname of the bai
  """
  # get ls of all files folder
  samples = os.popen('gsutil -m ls -al ' + gsfolder + '**.bai').read().split('\n')
  # compute size filepath

  sizes = {'gs://' + val.split('gs://')[1].split('#')[0]: int(val.split("2019-")[0]) for val in samples[:-2]}
  names = {}
  for k, val in sizes.items():
    if val in names:
      names[val].append(k)
    else:
      names[val] = [k]
  # get all bai in tsv
  samp = dm.WorkspaceManager(workspace).get_samples()
  for k, val in samp.iterrows():
    if val[bamcol] != 'NA' and val[baicol] != 'NA':
      # if bai has duplicate size
      code = os.system('gsutil ls ' + val[bamcol])
      if code == 256:
        if val[bamcol] is None:
          print('we dont have bam value for ' + str(k))
          continue
        else:
          print('no match values for ' + str(val[bamcol]))

        for va in names[sizes[val[baicol]]]:
          # for all duplicate size
          # if ls bam of bai duplicate size work
          # mv bam to bampath in gsfolder
          if '.bam' in va:
            if os.system('gsutil ls ' + va.split('.bam.bai')[0] + ".bam") == 0:
              print('gsutil mv ' + va.split('.bam.bai')[0] + ".bam " + val[bamcol])
              os.system('gsutil mv ' + va.split('.bam.bai')[0] + ".bam " + val[bamcol])
              break
          elif os.system('gsutil ls ' + va.split('.bai')[0] + ".bam") == 0:
            print('gsutil mv ' + va.split('.bai')[0] + ".bam " + val[bamcol])
            os.system('gsutil mv ' + va.split('.bai')[0] + ".bam " + val[bamcol])
            break
      elif code == signal.SIGINT:
        print('Awakened')
        break
    else:
      print("no data for " + str(k))


def shareTerraBams(users, workspace, samples, bamcols=["internal_bam_filepath", "internal_bai_filepath"]):
  """
  will share some files from gcp with a set of users using terra as metadata repo.

  only works with files that are listed on a terra workspace tsv but actually
  point to a regular google bucket and not a terra bucket.

  Args:
  ----
    users: list[str] of users' google accounts
    workspace: str namespace/workspace from url typically
    samples list[str] of samples_id for which you want to share data
    bamcols: list[str] list of column names of gsfiles to share

  Returns:
  --------
    a list of the gs path we have been giving access to
  """
  if type(users) is str:
    users = [users]
  wm = dm.WorkspaceManager(workspace)
  togiveaccess = np.ravel(wm.get_samples()[bamcols].loc[samples].values)
  for user in users:
    files = ''
    for i in togiveaccess:
      files += ' ' + i
    code = os.system("gsutil -m acl ch -ru " + user + ":R " + files)
    if code == signal.SIGINT:
      print('Awakened')
      break
  print('the files are stored here:\n\n')
  print(togiveaccess)
  print('\n\njust install and use gsutil to copy them')
  print('https://cloud.google.com/storage/docs/gsutil_install')
  print('https://cloud.google.com/storage/docs/gsutil/commands/cp')
  return togiveaccess


def shareCCLEbams(users, samples, groups=[], raise_error=True, arg_max_length=100000, bamcols=["internal_bam_filepath", "internal_bai_filepath"],
                  refsheet_url="https://docs.google.com/spreadsheets/d/1XkZypRuOEXzNLxVk9EOHeWRE98Z8_DBvL4PovyM01FE",
                  privacy_sheeturl="https://docs.google.com/spreadsheets/d/115TUgA1t_mD32SnWAGpW9OKmJ2W5WYAOs3SuSdedpX4"):
  """
  same as shareTerraBams but is completed to work with CCLE bams from the CCLE sample tracker

  You need to have gsheet installed and you '~/.client_secret.json', '~/.storage.json' set up

  Args:
  ----
    users: list[str] of users' google accounts
    groups: list[str] of groups' google accounts
    samples list[str] of samples cds_ids for which you want to share data
    bamcols: list[str] list of column names where bams/bais are
    raise_error: whether or not to raise an error if we find blacklisted lines
    refsheet_url: the google spreadsheet where the samples are stored
    privacy_sheeturl: the google spreadsheet where the samples are stored

  Returns:
  --------
    a list of the gs path we have been giving access to
  """
  sheets = Sheets.from_files('~/.client_secret.json', '~/.storage.json')
  print("You need to have gsheet installed and you '~/.client_secret.json', '~/.storage.json' set up")
  privacy = sheets.get(privacy_sheeturl).sheets[6].to_frame()
  refdata = sheets.get(refsheet_url).sheets[0].to_frame(index_col=0)
  blacklist = [i for i in privacy['blacklist'].values.tolist() if i is not np.nan]
  blacklisted = set(blacklist) & set(samples)
  print("we have " + str(len(blacklist)) + ' blacklisted files')
  if len(blacklisted):
    print("these lines are blacklisted " + blacklisted)
    if raise_error:
      raise ValueError("blacklistedlines")
  if type(users) is str:
    users = [users]

  togiveaccess = np.ravel(refdata[bamcols].loc[samples].values)
  usrs = ""
  for group in groups:
    usrs += " -g " + group + ":R"
  for user in users:
    usrs += " -u " + user + ":R"
  cmd_prefix = "gsutil -m acl ch " + usrs
  cmd = cmd_prefix
  for n, filename in enumerate(togiveaccess):
    oldcmd = cmd
    cmd += ' ' + filename
    if (len(cmd) > arg_max_length) | (n==len(togiveaccess)-1):
        if n < len(togiveaccess)-1:
            cmd = oldcmd
        print('granting access to {:d} files'.format(n))
        with open('/tmp/grantaccess{:d}.sh'.format(n), 'w') as f:
          f.write(cmd)
        code = os.system(cmd)
        cmd = cmd_prefix + ' ' + filename
        if code == signal.SIGINT:
          print('Awakened')
          return

  print('the files are stored here:\n\n' + refsheet_url)
  print('\n\njust install and use gsutil to copy them')
  print('https://cloud.google.com/storage/docs/gsutil_install')
  print('https://cloud.google.com/storage/docs/gsutil/commands/cp')
  return togiveaccess


def saveConfigs(workspace, filepath):
  """
  will save everything about a workspace into a csv and json file

  Args:
  -----
    workspace: str namespace/workspace from url typically
      namespace (str): project to which workspace belongs
      workspace (str): Workspace name
    filepath to save files
  """
  wm = dm.WorkspaceManager(workspace)
  h.createFoldersFor(filepath)

  conf = wm.get_configs()
  conf.to_csv(filepath + '.csv')
  params = {}
  params['GENERAL'] = wm.get_workspace_metadata()
  for k, val in conf.iterrows():
    params[k] = wm.get_config(val['name'])
  h.dictToFile(params, filepath + '.json')


def cleanWorkspace(workspaceid, only=[], toleave=[], defaulttoleave=['workspace', 'scripts',
                                                                     'notebooks', 'files', 'data', 'hound', 'references', 'name', 'folder']):
  """
  removes all processing folder in a terra workspace easily

  args:
    workspaceid: str, the workspace
    toleave: a list of first order folder in the bucket that you don't want to be deleted
    defaulttoleave: it should contain non processing folders that contain metadata and files for the workspace
  """
  toleave.extend(defaulttoleave)
  bucket = dm.WorkspaceManager(workspaceid).get_bucket_id()
  res = subprocess.run('gsutil -m ls gs://' + bucket, shell=True, capture_output=True)
  if res.returncode != 0:
    raise ValueError(str(res.stderr))
  res = str(res.stdout)[2:-1].split('\\n')[:-1]
  toremove = [val for val in res if val.split('/')[-2] not in toleave]
  if only:  # you were here
    toremove = [val for val in res if val.split('/')[-2] in only]
  if h.askif('we are going to remove ' + str(len(toremove)) + " files/folders:\n" + str(toremove) + "\nare you sure?"):
    gcp.rmFiles(toremove, add='-r')
  else:
    print("aborting")


def changeToBucket(samples, gsfolderto, name_col=None, values=['bam', 'bai'], filetypes=None, catchdup=False,
                   test=True):
  """
  moves all bam/bai files in a sampleList from Terra, to another gs bucket and rename them in the sample list

  will prevent erasing a duplicate sample by adding a random string or by flagging them and not copying them

  Args:
  ----
    samples: pandas.dataframe with columns to move
    gsfolderto: the bucket path to move the data to
    values: list of the cols in the dataframe containing the gs object path to be moved
    filetypes: list[str] of size values for each columns, give a suffix (.txt, .bam, ...)
    catchdup: if false will prepend a random string to the names before moving them, else will flag duplicate names
    test: only shows the output but does not move the files

  Returns:
  --------
    the updated sample pandas.dataframe
  """
  # to do the download to the new dataspace
  for i, val in samples.iterrows():
    ran = h.randomString(6, 'underscore', withdigits=False)
    for j, ntype in enumerate(values):
      # TODO try:catch
      filetype = '.'.join(val[ntype].split(
          '/')[-1].split('.')[1:]) if filetypes is None else filetypes[j]
      if name_col is None:
        name = val[ntype[0]].split('/')[-1].split('.')[0]
      elif name_col == "index":
        name = val.name
      else:
        name = val[name_col]
      name = name + '.' + filetype if catchdup else name + '_' + ran + '.' + filetype
      if not gcp.exists(gsfolderto + name) or not catchdup:
        cmd = 'gsutil cp ' + val[ntype] + ' ' + gsfolderto + name
        if test:
          print(cmd)
        else:
          res = subprocess.run(cmd, shell=True, capture_output=True)
          if res.returncode != 0:
            raise ValueError(str(res.stderr))
          samples.loc[i, ntype] = gsfolderto + name
      else:
        print(name + ' already exists in the folder: ' + gsfolderto)
        print(gcp.lsFiles([gsfolderto + name], '-la'))
  return samples


def delete_job(workspaceid, subid, taskid, DeleteCurrent=False, dryrun=True):
    """
    removes files generated by a job on Terra

    Args:
    -----
      workspaceid: str wokspace name
      subid: str the name of the job
      taskid: str the name of the task in this job
      DeleteCurrent: bool whether or not to delete files if they appear in one of the sample/samplesets/pairs data tables
      dryrun: bool just plot the commands but don't execute them
    """
    wm = dm.WorkspaceManager(workspaceid)
    bucket = wm.get_bucket_id()
    data= []
    if DeleteCurrent:
      if dryrun:
          print('gsutil -m rm gs://'+bucket+'/'+subid+'/*/'+taskid+'/**')
      else:
          res = subprocess.run('gsutil -m rm gs://'+bucket+'/'+subid+'/*/'+taskid+'/**', shell=True, capture_output=True)
          if res.returncode != 0:
              raise ValueError(str(res.stderr))
    else:
      res = subprocess.run('gsutil -m ls gs://'+bucket+'/'+subid+'/*/'+taskid+'/**', shell=True, capture_output=True)
      if res.returncode != 0 or len(str(res.stdout)) < 4:
          raise ValueError(str(res.stderr))
      data += str(res.stdout)[2:-1].split('\\n')[:-1]
      if "TOTAL:" in data[-1]:
          data = data[:-1]
      sam = pd.concat([wm.get_samples(), wm.get_pairs(), wm.get_sample_sets()])
      tokeep = set([val for val in sam.values.ravel() if type(val) is str and val[:5]=='gs://'])
      torm = set(data) - tokeep
      if dryrun:
        print(torm)
      else:
        h.parrun(['gsutil rm '+i for i in torm], cores=12)


#removing things from old failed workflows
def removeFromFailedWorkflows(workspaceid, maxtime = '2020-06-10', everythingFor=[], dryrun=False):
    """
    Lists all files from all jobs that have failed and deletes them.

    Can be very long

    Args:
    -----
      workspaceid: str the workspace name
      maxtime: str date format (eg. 2020-06-10) does not delete files generated past this date
      everythingFor: list[str] removes from these workflows even if not failed
      dryrun: bool whether or not to execute or just print commands
    """
    wm = dm.WorkspaceManager(workspaceid)
    for k, val in wm.get_submission_status(filter_active=False).iterrows():
        if (val.Failed > 0 or val.configuration in everythingFor) and val.date.date() > pd.to_datetime(maxtime):
            for w in wm.get_submission(val.submission_id)['workflows']:
                if w['status']=='Failed' or val.configuration in everythingFor:
                    try:
                        a = w['workflowId']
                    #else it was not even run
                    except:
                        continue
                    delete_job(workspaceid,val.submission_id,a,dryrun=dryrun)


def deleteHeavyFiles(workspaceid, unusedOnly=True):
    """
    deletes all files above a certain size in a workspace (that are used or unused)

    Args:
    ----
      workspaceid: str the name off the workspace
      unusedOnly: bool whether to delete used files as well (files that appear in one of the sample/samplesets/pairs data tables)
    """
    wm = dm.WorkspaceManager(workspaceid)
    bucket = wm.get_bucket_id()
    sizes = gcp.get_all_sizes('gs://'+bucket+'/')
    print('we got '+str(len(sizes))+' files')
    a  = list(sizes.keys())
    a.sort()
    ma = 100
    torm = []
    tot = 0
    for i in a[::-1]:
      if i>1000000*ma:
        tot += i
        for val in sizes[i]:
          torm.append(val)
    print('we might remove more than '+str(tot/1000000000)+'GB')
    if unusedOnly:
      sam = pd.concat([wm.get_samples(),wm.get_pairs(),wm.get_sample_sets()])
      tokeep = set([val for val in sam.values.ravel() if type(val) is str and val[:5]=='gs://'])
      torm = set(torm) - tokeep
    return torm


def findFilesInWorkspaces(names=[], lookup=['**', '*.', '.*']):
    """
    given All your terra workspaces, find a given gs filename

    Args:
    -----
      names: list[str] of filenames to find
      lookup: list[str] a set of flags giving how to look
        [** through all folders, *. can be preprended with anything,
        .* can be appended with anything]
    """
    ws = dm.list_workspaces()
    print('listing workspacs')
    file = []
    res = []
    for val in ws:
        val = val['workspace']
        print(val['namespace']+"/"+val['name'])
        buck = 'gs://'+val['bucketName']+"/"
        if subprocess.run('gsutil ls '+buck, capture_output=True, shell=True).returncode != 0:
            print("cannot access this bucket")
            continue
        if len(names) == 0:
            file.append(buck)
        if '**' in lookup:
            buck += '**'
        if not '*.' in lookup:
            buck += '/'
        for name in names:
            val = buck+name
            if '.*' in lookup:
                val += '*'
            data = subprocess.run("gsutil -m ls " + val,
                                  capture_output=True, shell=True)
            if data.returncode != 0:
                if "One or more URLs matched no objects" not in str(data.stderr):
                    raise ValueError(
                        'issue with the command: ' + str(data.stderr))
            if len(str(data.stdout)) < 4:
                continue
            res += str(data.stdout)[2:-1].split('\\n')[:-1]
            if "TOTAL:" in res[-1]:
                res = res[:-1]
    return res

def updateWorkflows(workflowIDs, path):
    """
    will download the latest version of workflows (from a list of workflowIDs) to a folder paths
    """

def uploadWorkflows(workspaceID, workflows, path=None):
    """
    updates the workflows on Terra and upgrades the workflow values on our workspace

    Args:
    -----
        workflows: dict(workflowID,location) or list(workflowID) if path
        path: folder path where files with same name as workflows' name are stored
    """
    method_folder="src/"
    methods = ['']
    dm.update_method()
