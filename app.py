from flask import Flask, render_template, request

from PubMedSearch import get_article_count_per_year, make_plot

app = Flask(__name__)


@app.route('/', methods=["POST", "GET"])
def homepage():
    """
    Pagina waar twee zoektermen kunnen worden opgegeven, vervolgens
    wordt per 5 jaar het aantal hits in pubmed in een barplot gezet.
    :return:
    """
    template_kwargs = {}
    if request.method == "POST":
        template_kwargs |= {**request.form}

        search_results = {}
        if term1 := request.form["search_term1"]:
            search_results[term1] = get_article_count_per_year(term1)
        if term2 := request.form["search_term2"]:
            search_results[term2] = get_article_count_per_year(term2)

        if search_results:
            plot_data = make_plot(search_results)
            template_kwargs["plot_data"] = plot_data

    return render_template("homepage.html", **template_kwargs)


if __name__ == '__main__':
    app.run()
