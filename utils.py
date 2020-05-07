from urllib import parse


data_dir = './static/data'


def format_query_string(data):
    result = ""

    for item in data:
        result += parse.quote_plus(item) + '=' + parse.quote_plus(data[item]) + '&'

    result = result[:-1]
    return result
