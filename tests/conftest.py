import os,shutil

def pytest_sessionstart(session):
     location = os.path.abspath(os.path.dirname(__file__))
#     shutil.rmtree(os.path.join(location,"data/testdata/kmer"), ignore_errors=True)
#     shutil.rmtree(os.path.join(location,"service/testdata/kmer"), ignore_errors=True)
#     shutil.rmtree(os.path.join(location,"service/testdata/marker"), ignore_errors=True)
#     if os.path.isfile(os.path.join(location,"service/testdata/db.sqlite")):
#         os.remove(os.path.join(location,"service/testdata/db.sqlite"))
        
def pytest_collection_modifyitems(session, config, items):
    classOrder = ["GeneralTestCase","IndexTestCase","ServiceDataTestCase","ServiceTestCase"]    
    classMapping = {item: item.cls.__name__ for item in items}
    #sort
    sortedItems = items.copy()
    for class_ in classOrder:
        sortedItems = [it for it in sortedItems if classMapping[it] != class_] + [
            it for it in sortedItems if classMapping[it] == class_
        ]
    items[:] = sortedItems


