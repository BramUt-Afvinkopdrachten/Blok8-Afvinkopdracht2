import base64
import numpy as np
from io import BytesIO

from Bio import Entrez
from matplotlib.figure import Figure

Entrez.email = ""


def get_article_count(search_term: str, mindate: str, maxdate: str) -> str:
    """
    Zoekt in pubmed het aantal artikelen waarin een zoekterm voorkomt.
    :param search_term: Zoekterm
    :param mindate: Vroegste datum
    :param maxdate: Laatste datum.
    :return: Aantal hits.
    """
    handle = Entrez.esearch(db="pubmed",
                            term=search_term,
                            rettype="count",
                            mindate=mindate,
                            maxdate=maxdate)
    record = Entrez.read(handle)
    handle.close()
    if count := record["Count"]:
        return count
    raise Exception("Something went wrong.")


def get_article_count_per_year(search_term: str):
    """
    :param str search_term: Zoekterm
    :return: List met tuples met per periode het aantal hits.
    vb. [('1971-1976', 0), ...]
    """
    minyear = 1971
    maxyear = 2021
    data = []
    for year in range(minyear, maxyear + 1, 5):
        count = get_article_count(search_term, str(year), str(year + 5))
        data.append((f"{year}-{year+5}", int(count)))
    return data


def make_plot(search_results: dict[str, tuple]):
    """
    Maakt een bar plot voor het aantaal hits op pubmed per periode.
    :param dict[str, tuple] search_results: dictionary met de zoektermen
    als keys en de jaren met aantal hits in een tuple als value.
    vb. search_results = {'covid': [('1971-1976', 0), ...]}
    :return: Image data van de plot.
    """
    fig = Figure()
    ax = fig.subplots()
    for c, (term, results) in enumerate(search_results.items()):
        years, counts = zip(*results)
        x = np.arange(len(years))
        ax.bar(x + (c * 0.25), counts, width=0.25, label=term)
    ax.set_xticks(x)
    ax.set_xticklabels(years)
    fig.autofmt_xdate()

    ax.legend()
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    data = base64.b64encode(buf.getbuffer()).decode("ascii")
    return data


if __name__ == '__main__':
    print(get_article_count_per_year("toxine"))
