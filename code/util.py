def prettify_medication_name(name, length=50):
  tokens = name.split("_")
  abbr_route = {
    "inhalation":"(Inhal)",
    "oral":"(Oral)",
    "intravenous":"(IV)",
    "nasal":"(Nas)",
    "injection":"(Inj)",
    "topical":"(Top)",
    "ophthalmic":"(Oph)",
    "subcutaneous":"(Subc)",
    "rectal":"(Rec)",
    "misc.(non-drug; combo route)":""
  }
  return abbr(tokens[0], length).title() + ((" "+abbr_route[tokens[1].lower()]) if tokens[1]!="unknown" else "")

def prettify_compound_medication_name(name, length=50, mode="all"):
  if mode=="trunc":
    len_first_name = len(prettify_medication_name(name.split("|")[0], length))
    return abbr("; ".join([prettify_medication_name(x, length) for x in name.split("|")]), max(length,len_first_name+2))
  if mode=="share":
    names = name.split("|")
    lengths = [int(np.ceil(length * (len(n) / len(name)))) for n in names]
    return "; ".join([prettify_medication_name(n, length) for n,length in zip(names,lengths)])
  if mode=="all":
    return "; ".join([prettify_medication_name(x, length) for x in name.split("|")])

def abbr(text, n=25):
  if len(text)>n-2:
    return text[:n-2]+'...'
  else:
    return text

class ProgressBar:
  def __init__(self, length):
    self.text =''
    self.length = length
    self.start_time = time.time()
  def update(self, progress, caption=''):
    self.progress = progress
    print('\b'*len(self.text), end='', flush=True)
    elapsed = time.time() - self.start_time
    estimated = elapsed * (self.length / max(1, self.progress))
    self.text = '%.2fs/%.2fs: (%.1f/%.1f) %s'%(elapsed, estimated, self.progress, self.length, caption)
    print(self.text, end='', flush=True)
